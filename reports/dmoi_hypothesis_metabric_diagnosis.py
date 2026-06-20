#!/usr/bin/env python3
"""Why the ER-lineage hypothesis is SUPPORTED in TCGA but REFUTED in METABRIC — a root-cause diagnosis.

The hypothesis screen measures a CONDITIONAL association ("does the hypothesis add beyond the textbook
anchor?"), i.e. the partial correlation of the hypothesis with the endpoint controlling for the anchor. That
quantity is sensitive to the joint correlation structure of anchor and hypothesis. This script computes the
same diagnostics in both cohorts to localize the cause of the non-reproduction.

Finding: estrogen-response separates LumA from LumB in BOTH cohorts (Cohen d ~ -0.2). The difference is the
proliferation-ER correlation: +0.19 in TCGA vs -0.17 in METABRIC (a sign flip). In METABRIC, LumB's lower ER
is collinear with its higher proliferation, so adjusting for proliferation removes the ER signal (partial
corr ~ 0); in TCGA the two are not collinear, so ER carries orthogonal information (partial corr -0.24). Hence
ER is SUPPORTED in TCGA and REFUTED in METABRIC -- not because ER biology is absent, but because its
information is redundant with proliferation in METABRIC's data structure. Writes
hypothesis_metabric_diagnosis.csv.

Run:  BRCA_DIR=/path/to/tcga_brca METABRIC_DIR=/path/to/metabric python reports/dmoi_hypothesis_metabric_diagnosis.py
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo
from sklearn.metrics import roc_auc_score
B = os.environ.get("BRCA_DIR", ""); M = os.environ.get("METABRIC_DIR", "")
PROLIF = ["MKI67", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB", "BUB1", "CCNE1", "CDC20",
          "TOP2A", "TYMS", "RRM2", "UBE2C", "CENPF", "FOXM1", "MELK", "KIF2C", "NUSAP1", "PTTG1"]
ER = ["ESR1", "GATA3", "FOXA1", "XBP1", "TFF1", "PGR", "GREB1", "CA12", "SLC39A6", "NAT1", "AR", "MLPH"]


def _cohen_d(x, g):
    a, b = x[g == 0], x[g == 1]
    sp = np.sqrt((a.std(ddof=1) ** 2 + b.std(ddof=1) ** 2) / 2) + 1e-9
    return (b.mean() - a.mean()) / sp


def _diag(name, E, y):
    P = mo.signature_score(E.loc[[g for g in PROLIF if g in E.index]].T, [g for g in PROLIF if g in E.index])
    R = mo.signature_score(E.loc[[g for g in ER if g in E.index]].T, [g for g in ER if g in E.index])
    Pz = (P - P.mean()) / P.std(); Rz = (R - R.mean()) / R.std(); yz = y.astype(float) - y.mean()

    def resid(v, ctrl):
        b = np.polyfit(ctrl, v, 1); return v - (b[0] * ctrl + b[1])
    pc = float(np.corrcoef(resid(Rz, Pz), resid(yz, Pz))[0, 1])
    ht = mo.hypothesis_anchor_test(P, R, y, cv=4, inner_repeats=1)   # commonality + mediation re-characterization
    return dict(cohort=name, n=int(len(y)), lumB=int(y.sum()),
               prolif_auroc=round(max(roc_auc_score(y, P), 1 - roc_auc_score(y, P)), 3),
               er_auroc=round(max(roc_auc_score(y, R), 1 - roc_auc_score(y, R)), 3),
               corr_prolif_er=round(float(np.corrcoef(Pz, Rz)[0, 1]), 3),
               partial_corr_er_given_prolif=round(pc, 3),
               cohen_d_prolif=round(float(_cohen_d(P, y)), 2), cohen_d_er=round(float(_cohen_d(R, y)), 2),
               unique_er_r2=ht["unique_hypothesis_r2"], common_r2=ht["common_r2"], redundancy=ht["redundancy"],
               prop_mediated=ht["prop_mediated"], collinearity_label=ht["collinearity_label"])


def main():
    assert B and os.path.isdir(B) and M and os.path.isdir(M), "set BRCA_DIR and METABRIC_DIR"
    ex = pd.read_csv(os.path.join(B, "HiSeqV2.gz"), sep="\t", index_col=0); ex = ex[~ex.index.duplicated()]
    cl = pd.read_csv(os.path.join(B, "BRCA_clinicalMatrix.tsv"), sep="\t", index_col=0)
    pam = cl["PAM50Call_RNAseq"].reindex(ex.columns); m = pam.isin(["LumA", "LumB"])
    r1 = _diag("TCGA", ex.loc[:, m], (pam[m] == "LumB").astype(int).values)
    mr = pd.read_csv(os.path.join(M, "mrna_microarray.txt"), sep="\t").drop(
        columns=["Entrez_Gene_Id"]).drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol")
    cp = pd.read_csv(os.path.join(M, "clinical_patient.txt"), sep="\t", comment="#")
    sub = cp.set_index(cp.columns[0])["CLAUDIN_SUBTYPE"]
    s = [c for c in mr.columns if c in sub.index and sub[c] in ("LumA", "LumB")]
    y = np.array([1 if sub[c] == "LumB" else 0 for c in s])
    Em = mr[s].apply(pd.to_numeric, errors="coerce"); Em = Em.T.fillna(Em.mean(1)).T
    r2 = _diag("METABRIC", Em, y)
    out = pd.DataFrame([r1, r2])
    out.to_csv(os.path.join(REPO, "hypothesis_metabric_diagnosis.csv"), index=False)
    print(out.to_string(index=False))
    print("\nROOT CAUSE: ER separates LumA/LumB in both cohorts (cohen_d_er ~ -0.2), but corr(prolif,ER) flips "
          "sign (+0.19 TCGA vs -0.17 METABRIC). In METABRIC ER is collinear with proliferation, so the "
          "partial correlation (ER given prolif) collapses to ~0 -> REFUTED; in TCGA it stays -0.24 -> SUPPORTED.")


if __name__ == "__main__":
    main()
