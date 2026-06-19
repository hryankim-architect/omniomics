#!/usr/bin/env python3
"""Clinical significance of the discovered basal axis: an honest separation of identity from outcome.

A reproducible biological axis need not be prognostic. We test whether the breast basal/keratinization
signature score is associated with overall survival in TCGA-BRCA, both alone and adjusted for the
proliferation anchor (the discovery's reference), and whether it tracks the clinical group it should
(ER-negative / basal-like disease).

Result: the basal score is NOT associated with overall survival (univariate Cox HR ~1.01, p~0.89; adjusted
for proliferation p~0.36; KM median-split log-rank p~0.29), while proliferation IS (p~0.001). But the basal
score does mark ER-negative/basal-like disease (AUROC ~0.70). So the discovered axis captures lineage
IDENTITY, not outcome -- reported honestly. Writes clinical_basal_survival.csv.

Run:  BRCA_DIR=/path/to/tcga_brca python reports/dmoi_clinical_survival.py
Deps: lifelines, scikit-learn.
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from sklearn.metrics import roc_auc_score
B = os.environ.get("BRCA_DIR", "")
PROLIF = ["MKI67", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB", "BUB1", "CCNE1", "CDC20",
          "TOP2A", "TYMS", "RRM2", "UBE2C", "CENPF", "FOXM1", "MELK", "KIF2C", "NUSAP1", "PTTG1"]


def main():
    assert B and os.path.isdir(B), "set BRCA_DIR to the tcga_brca folder (HiSeqV2.gz + BRCA_clinicalMatrix.tsv)"
    expr = pd.read_csv(os.path.join(B, "HiSeqV2.gz"), sep="\t", index_col=0); expr = expr[~expr.index.duplicated()]
    cl = pd.read_csv(os.path.join(B, "BRCA_clinicalMatrix.tsv"), sep="\t", index_col=0)
    basal = list(pd.read_csv(os.path.join(REPO, "novel_genes.csv"))["gene"])

    def score(genes):
        m = expr.loc[[g for g in genes if g in expr.index]]
        return m.sub(m.mean(1), axis=0).div(m.std(1) + 1e-9, axis=0).mean(0)

    bs, ps = score(basal), score(PROLIF)
    ot = pd.to_numeric(cl["OS_Time_nature2012"], errors="coerce"); oe = pd.to_numeric(cl["OS_event_nature2012"], errors="coerce")
    df = pd.DataFrame({"basal": bs, "prolif": ps}).join(pd.DataFrame({"T": ot, "E": oe})).dropna()
    df = df[df["T"] > 0]
    for c in ["basal", "prolif"]:
        df[c] = (df[c] - df[c].mean()) / df[c].std()
    uni = CoxPHFitter().fit(df[["basal", "T", "E"]], "T", "E").summary.loc["basal"]
    adj = CoxPHFitter().fit(df[["basal", "prolif", "T", "E"]], "T", "E").summary
    hb, hp = adj.loc["basal"], adj.loc["prolif"]
    hi = df["basal"] > df["basal"].median()
    lr = logrank_test(df["T"][hi], df["T"][~hi], df["E"][hi], df["E"][~hi])
    ercol = [c for c in cl.columns if c.lower() in ("er_status_by_ihc", "breast_carcinoma_estrogen_receptor_status")]
    er = (cl[ercol[0]] if ercol else cl.get("ER_Status_nature2012")).reindex(bs.index)
    mask = er.isin(["Negative", "Positive"]); y = (er[mask] == "Negative").astype(int).values
    auc_er = roc_auc_score(y, bs[mask].values)
    pd.DataFrame([
        dict(cohort="TCGA-BRCA", endpoint="overall_survival", n=int(len(df)), events=int(df["E"].sum()),
             basal_HR=round(float(np.exp(uni["coef"])), 3), basal_HR_p=round(float(uni["p"]), 4),
             basal_HR_adj_prolif=round(float(np.exp(hb["coef"])), 3), basal_adj_p=round(float(hb["p"]), 4),
             prolif_HR_adj=round(float(np.exp(hp["coef"])), 3), prolif_adj_p=round(float(hp["p"]), 4),
             km_logrank_p=round(float(lr.p_value), 4),
             note="basal score vs OS; adjusted model tests independence from the proliferation anchor"),
        dict(cohort="TCGA-BRCA", endpoint="ER_negative_status", n=int(mask.sum()), events=int(y.sum()),
             basal_HR="", basal_HR_p="", basal_HR_adj_prolif="", basal_adj_p="", prolif_HR_adj="", prolif_adj_p="",
             km_logrank_p="", note=f"clinical correlate (POSITIVE): basal score marks ER-negative/basal-like disease, AUROC={auc_er:.3f}"),
    ]).to_csv(os.path.join(REPO, "clinical_basal_survival.csv"), index=False)
    print(f"TCGA-BRCA OS n={len(df)} events={int(df['E'].sum())}")
    print(f"basal univariate HR={np.exp(uni['coef']):.2f} p={uni['p']:.3f} | adj p={hb['p']:.3f} | KM p={lr.p_value:.3f}")
    print(f"proliferation (adj) HR={np.exp(hp['coef']):.2f} p={hp['p']:.3f}")
    print(f"basal vs ER-negative AUROC={auc_er:.3f}")
    print("CONCLUSION: the basal axis is a lineage/identity marker (ER-negative/basal-like), NOT an "
          "independent predictor of overall survival -- an honest negative.")


if __name__ == "__main__":
    main()
