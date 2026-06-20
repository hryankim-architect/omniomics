#!/usr/bin/env python3
"""Robust hypothesis screen: anchor-family averaging + multiple-testing control.

Strengthens reports/dmoi_hypothesis_screen.py two ways:
  1. ANCHOR-FAMILY AVERAGING — score each Hallmark hypothesis's signal beyond TWO proliferation anchors
     (the curated 20-gene index and the data-driven meta-PCNA), average the delta, and require agreement
     (both anchors must call it supported). This makes the ranking robust to the choice of textbook anchor.
  2. MULTIPLE-TESTING CONTROL — across the 50 hallmarks, convert mean delta to an empirical-null one-sided
     p-value (null spread from the redundant/negative hallmarks, MAD-based; no costly permutation) and apply
     Benjamini-Hochberg FDR. A hypothesis is robust-SUPPORTED only if mean delta exceeds the threshold under
     BOTH anchors AND survives FDR.

Writes hypothesis_screen_robust.csv. Run:
  BRCA_DIR=/path/to/tcga_brca MSIGDB_GMT=/path/to/h.all...symbols.gmt python reports/dmoi_hypothesis_screen_robust.py
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo
B = os.environ.get("BRCA_DIR", ""); GMT = os.environ.get("MSIGDB_GMT", "")
PROLIF = ["MKI67", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB", "BUB1", "CCNE1", "CDC20",
          "TOP2A", "TYMS", "RRM2", "UBE2C", "CENPF", "FOXM1", "MELK", "KIF2C", "NUSAP1", "PTTG1"]
ADD = 0.01; FDR_Q = 0.10


def main():
    assert B and os.path.isdir(B) and GMT and os.path.exists(GMT), "set BRCA_DIR and MSIGDB_GMT"
    expr = pd.read_csv(os.path.join(B, "HiSeqV2.gz"), sep="\t", index_col=0); expr = expr[~expr.index.duplicated()]
    cl = pd.read_csv(os.path.join(B, "BRCA_clinicalMatrix.tsv"), sep="\t", index_col=0)
    pam = cl["PAM50Call_RNAseq"].reindex(expr.columns); mask = pam.isin(["LumA", "LumB"])
    y = (pam[mask] == "LumB").astype(int).values; E = expr.loc[:, mask]
    # two proliferation anchors (the family)
    A1 = mo.signature_score(E.loc[[g for g in PROLIF if g in E.index]].T, [g for g in PROLIF if g in E.index])
    mp = mo.marker_correlated_anchor(E.T, marker="PCNA", top_k=50, exclude_marker=True)
    A2 = mo.signature_score(E.loc[mp].T, mp)
    # hallmark hypotheses
    hyp = {}
    with open(GMT) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t"); g = [x for x in parts[2:] if x in E.index]
            if len(g) >= 10:
                hyp[parts[0]] = mo.signature_score(E.loc[g].T, g)

    def delta(anchor, h):
        return float(mo.anchored_integrate(anchor.reshape(-1, 1), h.reshape(-1, 1), y, cv=4,
                                           random_state=0, inner_repeats=1)["delta"])
    from sklearn.metrics import roc_auc_score
    names = list(hyp); d1 = np.array([delta(A1, hyp[n]) for n in names])
    d2 = np.array([delta(A2, hyp[n]) for n in names])
    md = (d1 + d2) / 2
    auc = np.array([max(roc_auc_score(y, hyp[n]), 1 - roc_auc_score(y, hyp[n])) for n in names])
    # empirical-null FDR: null spread from MAD of mean deltas (dominated by redundant hallmarks ~ 0)
    sd0 = 1.4826 * np.median(np.abs(md - np.median(md))) + 1e-9
    p = norm.sf(md / sd0)
    q = mo.benjamini_hochberg(p)
    agree = (d1 > ADD) & (d2 > ADD)
    robust = agree & (q < FDR_Q)
    verdict = np.where(robust, "SUPPORTED",
                       np.where(auc >= 0.6, "EXPLAINED_BY_TEXTBOOK", "REFUTED"))
    out = pd.DataFrame(dict(hypothesis=names, delta_curated=np.round(d1, 4), delta_metaPCNA=np.round(d2, 4),
                            mean_delta=np.round(md, 4), both_anchors_support=agree, auroc_hypothesis=np.round(auc, 3),
                            q_value=np.round(q, 4), robust_verdict=verdict)).sort_values("mean_delta", ascending=False)
    out.to_csv(os.path.join(REPO, "hypothesis_screen_robust.csv"), index=False)
    print(f"robust screen: {int(robust.sum())} SUPPORTED (both anchors + FDR<{FDR_Q}) of {len(names)} hallmarks | null sd0={sd0:.4f}")
    print(out.head(10).to_string(index=False))
    print("\nproliferation hallmarks (expected EXPLAINED):")
    print(out[out["hypothesis"].str.contains("E2F_TARGETS|G2M_CHECKPOINT|MYC_TARGETS")][
        ["hypothesis", "mean_delta", "both_anchors_support", "robust_verdict"]].to_string(index=False))


if __name__ == "__main__":
    main()
