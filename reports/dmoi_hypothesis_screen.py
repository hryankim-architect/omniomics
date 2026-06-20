#!/usr/bin/env python3
"""Hypothesis screen: rank a standard pathway library against the textbook anchor + real data.

Generalizes the single hypothesis test into a screen — score every MSigDB Hallmark gene set as a candidate
"hypothesis anchor" and rank by the signal it adds BEYOND the textbook proliferation anchor on LumA-vs-LumB.
This turns the anchored frame into a principled, library-wide question: *which textbook pathways carry a real
axis beyond proliferation?* Per Venet et al. (2011), proliferation-type hallmarks should be redundant
(EXPLAINED_BY_TEXTBOOK) while genuinely orthogonal lineage programs should be SUPPORTED.

Writes hypothesis_screen_results.csv (one row per hallmark, sorted by delta beyond the textbook).

Run:  BRCA_DIR=/path/to/tcga_brca MSIGDB_GMT=/path/to/h.all...symbols.gmt python reports/dmoi_hypothesis_screen.py
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo
B = os.environ.get("BRCA_DIR", ""); GMT = os.environ.get("MSIGDB_GMT", "")
PROLIF = ["MKI67", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB", "BUB1", "CCNE1", "CDC20",
          "TOP2A", "TYMS", "RRM2", "UBE2C", "CENPF", "FOXM1", "MELK", "KIF2C", "NUSAP1", "PTTG1"]


def main():
    assert B and os.path.isdir(B) and GMT and os.path.exists(GMT), "set BRCA_DIR and MSIGDB_GMT"
    expr = pd.read_csv(os.path.join(B, "HiSeqV2.gz"), sep="\t", index_col=0); expr = expr[~expr.index.duplicated()]
    cl = pd.read_csv(os.path.join(B, "BRCA_clinicalMatrix.tsv"), sep="\t", index_col=0)
    pam = cl["PAM50Call_RNAseq"].reindex(expr.columns); mask = pam.isin(["LumA", "LumB"])
    y = (pam[mask] == "LumB").astype(int).values; E = expr.loc[:, mask]
    T = mo.signature_score(E.loc[[g for g in PROLIF if g in E.index]].T, [g for g in PROLIF if g in E.index])
    # build a hypothesis score from every Hallmark gene set
    hyp = {}
    with open(GMT) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t"); name = parts[0]
            genes = [g for g in parts[2:] if g in E.index]
            if len(genes) >= 10:
                hyp[name] = mo.signature_score(E.loc[genes].T, genes)
    ranked = mo.rank_hypotheses(T, hyp, y, cv=4, random_state=0, inner_repeats=1)
    df = pd.DataFrame([{"hypothesis": r["hypothesis"], "auroc_hypothesis": r["auroc_hypothesis"],
                        "corr_with_textbook": r["corr_textbook_hypothesis"],
                        "delta_beyond_textbook": r["delta_hyp_given_textbook"], "verdict": r["verdict"]}
                       for r in ranked])
    df.to_csv(os.path.join(REPO, "hypothesis_screen_results.csv"), index=False)
    nsup = int((df["verdict"] == "SUPPORTED").sum())
    print(f"screened {len(df)} Hallmark hypotheses vs the proliferation textbook anchor (LumA/B) | SUPPORTED={nsup}")
    print("\nTop hypotheses adding beyond the textbook:")
    print(df.head(8).to_string(index=False))
    prolif_like = df[df["hypothesis"].str.contains("E2F_TARGETS|G2M_CHECKPOINT|MYC_TARGETS", regex=True)]
    print("\nProliferation-type hallmarks (expected redundant with the textbook anchor):")
    print(prolif_like[["hypothesis", "auroc_hypothesis", "delta_beyond_textbook", "verdict"]].to_string(index=False))


if __name__ == "__main__":
    main()
