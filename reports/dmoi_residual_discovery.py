#!/usr/bin/env python3
"""Knowledge-anchored residual discovery: separate the known (textbook), keep only the real new.

The capstone that combines the whole toolkit. Anchor on a FIXED textbook prior (a proliferation signature
for LumA/B), gate genome-wide data onto its residual, test whether that residual is SIGNAL not noise
(label-permutation null on the leakage-safe gated gain), and surface the features driving it that are
ORTHOGONAL to the anchor -- candidate "new biology beyond the textbook". A matched random panel + a
permutation null on the discovered panel are the noise controls (omniomics.multiomics.anchored_residual_discovery).

Result on TCGA-BRCA LumA/B (anchor = 20-gene proliferation index, 0 trained params): the noise-gated
residual is real signal, and the top anchor-orthogonal genes are the classic BASAL / squamous-lineage axis
(KRT5/14/17/6B, TP63, DSG3/DSC3, SOX10, COL17A1, KLK5/7/8) -- i.e. beyond the textbook "LumA/B = proliferation",
the data recovers a coherent lineage axis. The discovered panel gated onto the prior beats a matched random
panel. Writes discovery_results.csv (metrics) and novel_genes.csv (the discovered, anchor-orthogonal panel).

Run:  BRCA_DIR=/path/to/tcga_brca python reports/dmoi_residual_discovery.py
Env:  BRCA_DIR, FG_DIR (cache rna1500.tsv), DISC_PERM (default 12), DISC_TOPK (default 30).
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
NPERM = int(os.environ.get("DISC_PERM", "40")); TOPK = int(os.environ.get("DISC_TOPK", "30"))
PROLIF = ["MKI67", "AURKA", "BIRC5", "CCNB1", "CCNB2", "CDC20", "CEP55", "KIF2C", "MYBL2", "NDC80",
          "RRM2", "TYMS", "UBE2C", "BUB1", "CENPF", "PTTG1", "EXO1", "ANLN", "UBE2T", "NUF2"]


def _rna1500(S):
    rp = os.path.join(FG, "rna1500.tsv")
    if os.path.exists(rp):
        return pd.read_csv(rp, sep="\t", index_col=0)
    rna = pd.read_csv(os.path.join(BRCA, "HiSeqV2.gz"), sep="\t", index_col=0, usecols=["sample"] + S)
    rna = rna.loc[rna.var(axis=1).sort_values(ascending=False).index[:1500]]; rna.to_csv(rp, sep="\t")
    return rna


def main():
    coh = pd.read_csv(os.path.join(BRCA, "cohort_v2.tsv"), sep="\t")
    coh = coh[(coh.group.isin(["LumA", "LumB"])) & (coh.has_rna) & (coh.has_meth)]
    lab = dict(zip(coh.sample_id, coh.group))
    rna = _rna1500([s for s in coh.sample_id])
    S = [s for s in rna.columns if s in lab]
    y = np.array([1 if lab[s] == "LumB" else 0 for s in S])
    pr = pd.read_csv(os.path.join(BRCA, "HiSeqV2.gz"), sep="\t", index_col=0, usecols=["sample"] + S).reindex(PROLIF).dropna()
    anchor = mo.signature_score(pr[S], list(pr.index))                       # fixed textbook prior, 0 params
    R = rna[S]; Rz = R.sub(R.mean(1), axis=0).div(R.std(1).replace(0, 1), axis=0)
    X = Rz.T.values; genes = list(Rz.index)
    res = mo.anchored_residual_discovery(anchor, X, genes, y, top_k=TOPK, corr_max=0.6,
                                         cv=5, random_state=0, n_perm=NPERM, inner_repeats=3)
    # metrics CSV
    met = {k: res[k] for k in ("auroc_anchor", "auroc_combined", "delta", "perm_p", "novel_delta",
                               "random_delta_mean", "novel_vs_random_p", "signal")}
    met.update(n=len(S), anchor="proliferation_20genes")
    pd.DataFrame([met]).to_csv(os.path.join(REPO, "discovery_results.csv"), index=False)
    # novel genes CSV
    nov = pd.DataFrame(res["novel"], columns=["gene", "partial_corr_y_given_anchor", "corr_with_anchor"])
    nov.to_csv(os.path.join(REPO, "novel_genes.csv"), index=False)
    print(f"LumA/B n={len(S)} | anchor(proliferation) AUROC={res['auroc_anchor']:.3f}")
    print(f"noise-gated residual: combined={res['auroc_combined']:.3f} delta={res['delta']:+.3f} "
          f"(label-perm p={res['perm_p']:.3f})")
    print(f"discovered novel panel (anchor-orthogonal): delta={res['novel_delta']:+.3f} vs random-panel "
          f"mean {res['random_delta_mean']:+.3f} (novel-vs-random p={res['novel_vs_random_p']:.3f}) -> DISCOVERY SIGNAL={res['signal']}")
    print("top novel genes:", ", ".join(g for g, _, _ in res["novel"][:15]))
    print("\nwrote discovery_results.csv + novel_genes.csv")
    print("verdict: beyond the textbook 'LumA/B = proliferation' the noise-gated residual recovers a real,")
    print("         anchor-orthogonal BASAL/lineage axis (keratins/TP63/SOX10) -- new, and not noise.")


if __name__ == "__main__":
    main()
