#!/usr/bin/env python3
"""Hypothesis screen on OTHER endpoints — HER2 and ER — to show the frame generalizes (and stays honest).

Repeats the Hallmark hypothesis screen for two more TCGA-BRCA endpoints, each against its own textbook anchor:
  - HER2 status, anchored on the ERBB2 17q12 amplicon (INCOMPLETE in TCGA, AUROC ~0.75) -> expect some
    hallmarks to add a real axis beyond the amplicon (SUPPORTED).
  - ER status, anchored on a textbook ER/luminal signature (COMPLETE, AUROC ~0.94) -> a specificity check:
    when the textbook anchor already saturates the endpoint, the screen should surface few/no SUPPORTED hits.

For each endpoint we rank the 50 MSigDB Hallmark sets by signal added beyond the anchor (gate onto the anchor
residual), with an empirical-null Benjamini-Hochberg FDR. Writes hypothesis_screen_her2.csv and
hypothesis_screen_er.csv.

Run:  BRCA_DIR=... MSIGDB_GMT=... python reports/dmoi_hypothesis_screen_endpoints.py
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo
from sklearn.metrics import roc_auc_score
B = os.environ.get("BRCA_DIR", ""); GMT = os.environ.get("MSIGDB_GMT", "")
ADD = 0.01; FDR_Q = 0.10
ENDPOINTS = {
    "her2": dict(col="HER2_Final_Status_nature2012", pos="Positive", neg="Negative",
                 anchor=["ERBB2", "GRB7", "STARD3", "PGAP3", "TCAP", "PNMT", "PSMD3", "GSDMB", "ORMDL3"]),
    "er":   dict(col="breast_carcinoma_estrogen_receptor_status", pos="Positive", neg="Negative",
                 anchor=["ESR1", "GATA3", "FOXA1", "XBP1", "TFF1", "PGR", "GREB1", "CA12", "SLC39A6", "NAT1", "AR", "MLPH"]),
}


def _hallmarks(E):
    hyp = {}
    with open(GMT) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t"); g = [x for x in parts[2:] if x in E.index]
            if len(g) >= 10:
                hyp[parts[0]] = mo.signature_score(E.loc[g].T, g)
    return hyp


def main():
    assert B and os.path.isdir(B) and GMT and os.path.exists(GMT), "set BRCA_DIR and MSIGDB_GMT"
    expr = pd.read_csv(os.path.join(B, "HiSeqV2.gz"), sep="\t", index_col=0); expr = expr[~expr.index.duplicated()]
    cl = pd.read_csv(os.path.join(B, "BRCA_clinicalMatrix.tsv"), sep="\t", index_col=0)
    for ep, cfg in ENDPOINTS.items():
        lab = cl[cfg["col"]].reindex(expr.columns)
        mask = lab.isin([cfg["pos"], cfg["neg"]]); y = (lab[mask] == cfg["pos"]).astype(int).values
        E = expr.loc[:, mask.values]
        ag = [g for g in cfg["anchor"] if g in E.index]
        T = mo.signature_score(E.loc[ag].T, ag)
        aT = roc_auc_score(y, T); aT = max(aT, 1 - aT)
        hyp = _hallmarks(E); names = list(hyp)
        d = np.array([float(mo.anchored_integrate(T.reshape(-1, 1), hyp[n].reshape(-1, 1), y, cv=4,
                                                  random_state=0, inner_repeats=1)["delta"]) for n in names])
        auc = np.array([max(roc_auc_score(y, hyp[n]), 1 - roc_auc_score(y, hyp[n])) for n in names])
        sd0 = 1.4826 * np.median(np.abs(d - np.median(d))) + 1e-9
        q = mo.benjamini_hochberg(norm.sf(d / sd0))
        verdict = np.where((d > ADD) & (q < FDR_Q), "SUPPORTED",
                           np.where(auc >= 0.6, "EXPLAINED_BY_TEXTBOOK", "REFUTED"))
        out = pd.DataFrame(dict(hypothesis=names, delta_beyond_textbook=np.round(d, 4),
                                auroc_hypothesis=np.round(auc, 3), q_value=np.round(q, 4), verdict=verdict)
                           ).sort_values("delta_beyond_textbook", ascending=False)
        out.to_csv(os.path.join(REPO, f"hypothesis_screen_{ep}.csv"), index=False)
        nsup = int((out["verdict"] == "SUPPORTED").sum())
        print(f"[{ep.upper()}] textbook anchor AUROC={aT:.3f} ({'incomplete' if aT < 0.85 else 'complete'}) | "
              f"SUPPORTED={nsup}/{len(names)}")
        print(out.head(5).to_string(index=False)); print()
    print("READING: an INCOMPLETE textbook anchor (HER2 amplicon) leaves room for SUPPORTED hits; a COMPLETE "
          "anchor (ER signature) should yield few/none — a specificity check on the screen.")


if __name__ == "__main__":
    main()
