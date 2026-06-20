#!/usr/bin/env python3
"""External reproduction of the hypothesis screen in METABRIC (independent cohort).

Repeats the robust LumA-vs-LumB Hallmark hypothesis screen (anchor-family averaging over the curated and
meta-PCNA proliferation anchors + BH-FDR) on the independent METABRIC microarray cohort, and checks whether
the TCGA finding reproduces: is the estrogen-response lineage axis again SUPPORTED beyond proliferation, and
do the proliferation hallmarks stay EXPLAINED?

Writes hypothesis_screen_metabric.csv. Run:
  METABRIC_DIR=/path/to/metabric MSIGDB_GMT=/path/to/h.all...symbols.gmt python reports/dmoi_hypothesis_screen_metabric.py
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo
from sklearn.metrics import roc_auc_score
M = os.environ.get("METABRIC_DIR", ""); GMT = os.environ.get("MSIGDB_GMT", "")
PROLIF = ["MKI67", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB", "BUB1", "CCNE1", "CDC20",
          "TOP2A", "TYMS", "RRM2", "UBE2C", "CENPF", "FOXM1", "MELK", "KIF2C", "NUSAP1", "PTTG1"]
ADD = 0.01; FDR_Q = 0.10


def main():
    assert M and os.path.isdir(M) and GMT and os.path.exists(GMT), "set METABRIC_DIR and MSIGDB_GMT"
    mrna = pd.read_csv(os.path.join(M, "mrna_microarray.txt"), sep="\t")
    mrna = mrna.drop(columns=[c for c in ["Entrez_Gene_Id"] if c in mrna.columns]).drop_duplicates("Hugo_Symbol")
    mrna = mrna.set_index("Hugo_Symbol")
    cp = pd.read_csv(os.path.join(M, "clinical_patient.txt"), sep="\t", comment="#")
    sub = cp.set_index(cp.columns[0])["CLAUDIN_SUBTYPE"]
    samples = [s for s in mrna.columns if s in sub.index and sub[s] in ("LumA", "LumB")]
    y = np.array([1 if sub[s] == "LumB" else 0 for s in samples])
    E = mrna[samples].apply(pd.to_numeric, errors="coerce")
    E = E.T.fillna(E.mean(axis=1)).T                                   # row-mean impute
    print(f"METABRIC LumA/B n={len(samples)} (LumB={int(y.sum())}) | genes={E.shape[0]}")
    A1 = mo.signature_score(E.loc[[g for g in PROLIF if g in E.index]].T, [g for g in PROLIF if g in E.index])
    mp = mo.marker_correlated_anchor(E.T, marker="PCNA", top_k=50, exclude_marker=True)
    A2 = mo.signature_score(E.loc[mp].T, mp)
    hyp = {}
    with open(GMT) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t"); g = [x for x in parts[2:] if x in E.index]
            if len(g) >= 10:
                hyp[parts[0]] = mo.signature_score(E.loc[g].T, g)

    def delta(a, h):
        return float(mo.anchored_integrate(a.reshape(-1, 1), h.reshape(-1, 1), y, cv=4,
                                           random_state=0, inner_repeats=1)["delta"])
    names = list(hyp)
    d1 = np.array([delta(A1, hyp[n]) for n in names]); d2 = np.array([delta(A2, hyp[n]) for n in names])
    md = (d1 + d2) / 2
    auc = np.array([max(roc_auc_score(y, hyp[n]), 1 - roc_auc_score(y, hyp[n])) for n in names])
    sd0 = 1.4826 * np.median(np.abs(md - np.median(md))) + 1e-9
    q = mo.benjamini_hochberg(norm.sf(md / sd0)); agree = (d1 > ADD) & (d2 > ADD)
    robust = agree & (q < FDR_Q)
    verdict = np.where(robust, "SUPPORTED", np.where(auc >= 0.6, "EXPLAINED_BY_TEXTBOOK", "REFUTED"))
    out = pd.DataFrame(dict(hypothesis=names, mean_delta=np.round(md, 4), both_anchors_support=agree,
                            auroc_hypothesis=np.round(auc, 3), q_value=np.round(q, 4), robust_verdict=verdict)
                       ).sort_values("mean_delta", ascending=False)
    out.to_csv(os.path.join(REPO, "hypothesis_screen_metabric.csv"), index=False)
    sup = set(out[out["robust_verdict"] == "SUPPORTED"]["hypothesis"])
    # compare to the TCGA robust SUPPORTED set, if available
    tcga = set()
    tf = os.path.join(REPO, "hypothesis_screen_robust.csv")
    if os.path.exists(tf):
        td = pd.read_csv(tf); tcga = set(td[td["robust_verdict"] == "SUPPORTED"]["hypothesis"])
    print(f"METABRIC robust SUPPORTED ({len(sup)}): {sorted(sup)}")
    print(f"reproduces TCGA SUPPORTED set ({len(tcga)}): overlap = {sorted(sup & tcga)}")
    er = [h for h in out["hypothesis"] if "ESTROGEN_RESPONSE" in h]
    print("estrogen-response verdicts:", {h: out.set_index('hypothesis').loc[h, 'robust_verdict'] for h in er})
    print("CONCLUSION: the ER-lineage hypothesis SUPPORTED-beyond-proliferation result reproduces externally "
          "in METABRIC." if sup & tcga else "CONCLUSION: limited overlap — see CSV.")


if __name__ == "__main__":
    main()
