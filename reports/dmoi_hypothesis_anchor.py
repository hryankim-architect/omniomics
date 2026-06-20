#!/usr/bin/env python3
"""Hypothesis-as-anchor: confirm / explain-away / refute a hypothesis against a textbook anchor + real data.

Extends the anchored frame from "robustness across equivalent textbook anchors" to "test a *non-equivalent*
hypothesis." A hypothesis is expressed as a candidate anchor (a gene set the researcher believes drives the
endpoint); `hypothesis_anchor_test` gates it onto the textbook anchor's residual and returns a 3-way verdict.
Per Venet et al. (2011): a hypothesis is a novel mechanism only if it ADDS beyond the dominant textbook prior
— not merely predicts on its own.

Demo on TCGA-BRCA LumA-vs-LumB with the textbook = 20-gene proliferation index:
  - basal/keratinization  (positive control — the known orthogonal axis)   -> expect SUPPORTED
  - immune / cytotoxic    (a plausible secondary hypothesis)               -> data decides
  - random gene set       (Venet-style negative control; seeded)           -> expect EXPLAINED_BY_TEXTBOOK / REFUTED

Writes hypothesis_anchor_results.csv.

Run:  BRCA_DIR=/path/to/tcga_brca python reports/dmoi_hypothesis_anchor.py
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo
B = os.environ.get("BRCA_DIR", "")
PROLIF = ["MKI67", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB", "BUB1", "CCNE1", "CDC20",
          "TOP2A", "TYMS", "RRM2", "UBE2C", "CENPF", "FOXM1", "MELK", "KIF2C", "NUSAP1", "PTTG1"]
HYP = {
    "basal_keratinization": ["KRT5", "KRT14", "KRT17", "KRT6B", "TP63", "DSG3", "DSC3", "SOX10", "COL17A1"],
    "immune_cytotoxic": ["CD8A", "CD3D", "GZMB", "PRF1", "IFNG", "CXCL9", "CXCL10", "GZMK", "NKG7", "CCL5"],
}


def main():
    assert B and os.path.isdir(B), "set BRCA_DIR"
    expr = pd.read_csv(os.path.join(B, "HiSeqV2.gz"), sep="\t", index_col=0); expr = expr[~expr.index.duplicated()]
    cl = pd.read_csv(os.path.join(B, "BRCA_clinicalMatrix.tsv"), sep="\t", index_col=0)
    pam = cl["PAM50Call_RNAseq"].reindex(expr.columns); mask = pam.isin(["LumA", "LumB"])
    y = (pam[mask] == "LumB").astype(int).values
    E = expr.loc[:, mask]
    T = mo.signature_score(E.loc[[g for g in PROLIF if g in E.index]].T, [g for g in PROLIF if g in E.index])
    # seeded random hypothesis (Venet-style negative control)
    rng = np.random.default_rng(0)
    HYP["random_30genes"] = list(rng.choice([g for g in E.index], 30, replace=False))
    rows = []
    for name, genes in HYP.items():
        g_in = [g for g in genes if g in E.index]
        Hh = mo.signature_score(E.loc[g_in].T, g_in)
        r = mo.hypothesis_anchor_test(T, Hh, y, cv=5, random_state=0, inner_repeats=2)
        rows.append(dict(hypothesis=name, n_genes=len(g_in), textbook="proliferation_20gene", **r))
        print(f"{name:22} aucH={r['auroc_hypothesis']:.3f} corr(T,H)={r['corr_textbook_hypothesis']:+.2f} "
              f"deltaH|T={r['delta_hyp_given_textbook']:+.4f} -> {r['verdict']}")
    pd.DataFrame(rows).to_csv(os.path.join(REPO, "hypothesis_anchor_results.csv"), index=False)
    print("\nwrote hypothesis_anchor_results.csv")
    print("READING: SUPPORTED = adds a real axis beyond the textbook; EXPLAINED_BY_TEXTBOOK = predicts alone "
          "but redundant once proliferation is controlled (Venet); REFUTED = neither.")


if __name__ == "__main__":
    main()
