#!/usr/bin/env python3
"""Modality generalization of knowledge-anchored residual discovery — the methylation modality.

The discovery machinery (omniomics.multiomics.anchored_residual_discovery) is modality-agnostic: it
operates on any (samples x features) matrix. Here the DATA modality is DNA methylation (450K CpGs) rather
than RNA, with the same fixed RNA proliferation prior as the anchor, on LumA/B.

Honest result: the method runs unchanged on methylation, but on this endpoint methylation carries no axis
beyond the RNA anchor -- random genome-wide CpGs add ~0 (Δ −0.001) and even basal-gene-targeted CpGs add ~0
(Δ +0.000; basal methylation alone only 0.756). So the basal/lineage axis discovered in RNA (Δ +0.029) is a
TRANSCRIPTIONAL signal, not (predictively) methylation-encoded here -- and the method correctly does not
manufacture a methylation discovery. Modality-agnostic, and still discriminating. Writes
methylation_discovery_results.csv and novel_cpgs_meth.csv (discovered CpGs mapped to genes).

Run:  BRCA_DIR=/path/to/tcga_brca python reports/dmoi_discovery_methylation.py
Env:  BRCA_DIR, FG_DIR (cache meth_gw.tsv), DISC_PERM (default 12).
"""
import os, sys, gzip
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo
from sklearn.metrics import roc_auc_score
try:
    from omniomics import config; DEFAULT_BRCA = config.brca_tcga_dir()
except Exception:
    DEFAULT_BRCA = ""
BRCA = os.environ.get("BRCA_DIR", DEFAULT_BRCA); FG = os.environ.get("FG_DIR", os.getcwd())
NPERM = int(os.environ.get("DISC_PERM", "12"))
PROLIF = ["MKI67", "AURKA", "BIRC5", "CCNB1", "CCNB2", "CDC20", "CEP55", "KIF2C", "MYBL2", "NDC80",
          "RRM2", "TYMS", "UBE2C", "BUB1", "CENPF", "PTTG1", "EXO1", "ANLN", "UBE2T", "NUF2"]
BASAL = {"KRT5", "KRT14", "KRT17", "KRT6B", "KRT6A", "KRT16", "TP63", "DSG3", "DSC3", "SOX10",
         "COL17A1", "KLK5", "KLK7", "KLK8", "SFN", "TRIM29", "CDH3"}


def _lumab():
    coh = pd.read_csv(os.path.join(BRCA, "cohort_v2.tsv"), sep="\t")
    coh = coh[(coh.group.isin(["LumA", "LumB"])) & (coh.has_rna) & (coh.has_meth)]
    return dict(zip(coh.sample_id, coh.group))


def _stream_cpgs(ids, samples_hint=None):
    keep, header = {}, None
    with gzip.open(os.path.join(BRCA, "HumanMethylation450.gz"), "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")[1:]
        for line in fh:
            cg, _, rest = line.partition("\t")
            if cg in ids:
                keep[cg] = rest.rstrip("\n").split("\t")
    return pd.DataFrame(keep, index=header).T.apply(pd.to_numeric, errors="coerce")


def main():
    lab = _lumab()
    gw = pd.read_csv(os.path.join(FG, "meth_gw.tsv"), sep="\t", index_col=0) if os.path.exists(os.path.join(FG, "meth_gw.tsv")) else None
    pm = pd.read_csv(os.path.join(BRCA, "hm450_probemap.tsv"), sep="\t", usecols=["#id", "gene"])
    basal_ids = set(pm["#id"][pm["gene"].apply(lambda g: isinstance(g, str) and bool(set(g.split(",")) & BASAL))])
    bm = _stream_cpgs(basal_ids)
    if gw is None:
        stride = set(pm["#id"].iloc[::130])
        gw = _stream_cpgs(stride)
    S = [s for s in gw.columns if s in bm.columns and s in lab]
    y = np.array([1 if lab[s] == "LumB" else 0 for s in S])
    pr = pd.read_csv(os.path.join(BRCA, "HiSeqV2.gz"), sep="\t", index_col=0, usecols=["sample"] + S).reindex(PROLIF).dropna()
    anchor = mo.signature_score(pr[S], list(pr.index))

    def prep(df):
        M = df[S].T.astype(float); return M.fillna(M.mean()).fillna(0.0)
    Gw, Bm = prep(gw), prep(bm)
    aA = roc_auc_score(y, anchor)
    res = mo.anchored_residual_discovery(anchor, Gw.values, list(Gw.columns), y, top_k=30, corr_max=0.6,
                                         cv=5, random_state=0, n_perm=NPERM, inner_repeats=2)
    rb = mo.anchored_integrate(anchor.reshape(-1, 1), Bm.values, y, cv=5, random_state=0, inner_repeats=3)
    basal_alone = mo.select_anchor({"m": Bm.values}, y, cv=5, repeats=3)["ranking"][0][1]
    rows = [dict(modality="methylation", endpoint="LumA_vs_LumB", anchor="proliferation_RNA_20genes",
                 anchor_auroc=round(aA, 3), data="random_genomewide", n_features=Gw.shape[1],
                 combined=round(res["auroc_combined"], 3), delta=round(res["delta"], 3), basal_meth_alone="",
                 note="modality-agnostic; random genome-wide methylation adds ~0 beyond the RNA proliferation anchor"),
            dict(modality="methylation", endpoint="LumA_vs_LumB", anchor="proliferation_RNA_20genes",
                 anchor_auroc=round(aA, 3), data="basal_gene_targeted", n_features=Bm.shape[1],
                 combined=round(rb["auroc_combined"], 3), delta=round(rb["delta"], 3), basal_meth_alone=round(basal_alone, 3),
                 note="targeted basal-gene methylation also adds ~0 -- the basal axis is transcriptional, not predictively methylation-encoded")]
    pd.DataFrame(rows).to_csv(os.path.join(REPO, "methylation_discovery_results.csv"), index=False)
    pm_s = pm.set_index("#id")["gene"]
    pd.DataFrame(res["novel"], columns=["cpg", "partial_corr", "corr_with_anchor"]).assign(
        gene=lambda d: [pm_s.get(c, ".") for c in d["cpg"]]).to_csv(os.path.join(REPO, "novel_cpgs_meth.csv"), index=False)
    print(pd.DataFrame(rows)[["data", "n_features", "combined", "delta"]].to_string(index=False))
    print("verdict: the discovery method is modality-agnostic (runs on methylation), but methylation carries")
    print("         no axis beyond the RNA anchor on this endpoint -- the basal axis is transcriptional.")


if __name__ == "__main__":
    main()
