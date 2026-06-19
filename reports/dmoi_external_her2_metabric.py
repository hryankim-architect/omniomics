#!/usr/bin/env python3
"""External validation of the HER2 discovery axis in METABRIC — an interpretable NEGATIVE.

Symmetric counterpart to reports/dmoi_external_metabric.py (which validated the basal axis). Reproduces the
HER2 amplicon-anchored residual discovery in METABRIC (microarray, HER2 IHC status, n=1980). Result: the
HER2 neuroendocrine/immune axis does NOT replicate -- because the ERBB2 amplicon anchor is near-complete in
METABRIC (AUROC ~0.997 vs 0.752 in TCGA: METABRIC HER2 calls are far more amplification-concordant), so
essentially no residual exists. The TCGA HER2 panel adds Δ ~0 (p=1.0) and an unbiased re-discovery overlaps
the TCGA panel only ~2/30 (vs 20/30 for basal). So reproducibility tracks biological coherence: the
pathway-enriched basal axis reproduces across cohorts, while the diffuse HER2 axis was a cohort-specific
residual; the method stays well-behaved (finds nothing where there is no residual, like TCGA ER). Writes
external_validation_her2_metabric.csv and novel_genes_her2_metabric.csv.

Run:  METABRIC_DIR=/path/to/metabric python reports/dmoi_external_her2_metabric.py
Env:  METABRIC_DIR, FG_DIR (cache mb_her2_expr.tsv / mb_her2_status.tsv).
"""
import os, sys
import numpy as np, pandas as pd
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo
from sklearn.metrics import roc_auc_score
MB = os.environ.get("METABRIC_DIR", ""); FG = os.environ.get("FG_DIR", os.getcwd())
AMP = ["ERBB2", "GRB7", "STARD3", "PGAP3", "TCAP", "PNMT", "PSMD3", "GSDMB", "ORMDL3"]


def _load():
    ep, lp = os.path.join(FG, "mb_her2_expr.tsv"), os.path.join(FG, "mb_her2_status.tsv")
    if os.path.exists(ep) and os.path.exists(lp):
        e = pd.read_csv(ep, sep="\t").drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol")
        if "Entrez_Gene_Id" in e.columns:
            e = e.drop(columns=["Entrez_Gene_Id"])
        return e, pd.read_csv(lp, sep="\t", index_col=0).iloc[:, 0]
    assert MB and os.path.isdir(MB), f"METABRIC dir not found: {MB!r} (set METABRIC_DIR)"
    cs = pd.read_csv(os.path.join(MB, "clinical_sample.txt"), sep="\t", comment="#")
    h = cs.set_index(cs.columns[0])["HER2_STATUS"]; h = h[h.isin(["Positive", "Negative"])]
    panel = list(pd.read_csv(os.path.join(REPO, "novel_genes_her2.csv"))["gene"])
    targets = set(AMP) | set(panel); keep = {}
    with open(os.path.join(MB, "mrna_microarray.txt")) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for i, line in enumerate(fh):
            g = line[:line.find("\t")]
            if g in targets or i % 25 == 0:
                keep[i] = line
    rows = [[header[0]] + header[2:]] + [[ln.rstrip("\n").split("\t")[0]] + ln.rstrip("\n").split("\t")[2:] for ln in keep.values()]
    e = pd.DataFrame(rows[1:], columns=rows[0]).drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol")
    return e.apply(pd.to_numeric, errors="coerce"), h


def main():
    e, h = _load()
    S = [s for s in e.columns if s in h.index]; y = np.array([1 if h[s] == "Positive" else 0 for s in S])
    e = e[S].apply(pd.to_numeric, errors="coerce")
    amp = [g for g in AMP if g in e.index]
    tcga = list(pd.read_csv(os.path.join(REPO, "novel_genes_her2.csv"))["gene"]); pan = [g for g in tcga if g in e.index]
    anchor = mo.signature_score(e.loc[amp].T, amp)

    def z(genes):
        M = e.loc[genes].T.astype(float); return M.fillna(M.mean()).fillna(0.0).values
    aA = roc_auc_score(y, anchor)
    rp = mo.anchored_integrate(anchor.reshape(-1, 1), z(pan), y, cv=5, random_state=0, inner_repeats=1)
    bg = [g for g in e.index if g not in set(amp) | set(pan)]; rng = np.random.default_rng(0)
    rand = np.array([mo.anchored_integrate(anchor.reshape(-1, 1), z(list(rng.choice(bg, len(pan), replace=False))),
                                           y, cv=5, random_state=0, inner_repeats=1)["delta"] for _ in range(15)])
    p_rep = (1 + int((rand >= rp["delta"]).sum())) / 16
    pool = [g for g in e.index if g not in set(amp)]
    res = mo.anchored_residual_discovery(anchor, z(pool), pool, y, top_k=30, corr_max=0.6, cv=5, random_state=0, n_perm=6, inner_repeats=1)
    overlap = [g for g, _, _ in res["novel"] if g in set(tcga)]
    pd.DataFrame(res["novel"], columns=["gene", "partial_corr", "corr_with_anchor"]).to_csv(
        os.path.join(REPO, "novel_genes_her2_metabric.csv"), index=False)
    out = pd.DataFrame([
        ["replication", "n", len(S), "METABRIC HER2 (independent expression cohort)"],
        ["replication", "amplicon_anchor_auroc_METABRIC", round(aA, 3), "near-complete here -- leaves ~no residual"],
        ["replication", "amplicon_anchor_auroc_TCGA", 0.752, "TCGA: amplicon incomplete -> a residual existed to discover"],
        ["replication", "tcga_her2_panel_delta", round(rp["delta"], 3), "TCGA HER2 panel adds ~0 in METABRIC"],
        ["replication", "random_panel_mean_delta", round(float(rand.mean()), 3), "matched random panels"],
        ["replication", "panel_vs_random_p", round(p_rep, 3), f"{int((rand>=rp['delta']).sum())}/15 random match"],
        ["rediscovery", "overlap_with_tcga_her2_of30", len(overlap), ", ".join(overlap) or "(none)"],
        ["verdict", "replicates", "No", "HER2 axis does NOT replicate -- amplicon near-complete in METABRIC (interpretable negative)"],
    ], columns=["test", "metric", "value", "note"])
    out.to_csv(os.path.join(REPO, "external_validation_her2_metabric.csv"), index=False)
    print(out.to_string(index=False))
    print(f"\nCONCLUSION: the HER2 axis does NOT externally validate -- in METABRIC the amplicon anchor is "
          f"near-complete (AUROC {aA:.3f}), so no residual exists. Reproducibility tracks coherence: basal "
          f"reproduces, diffuse HER2 is cohort-specific.")


if __name__ == "__main__":
    main()
