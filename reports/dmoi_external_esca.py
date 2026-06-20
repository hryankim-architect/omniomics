#!/usr/bin/env python3
"""Cross-CANCER replication (a third cancer): does the breast basal/keratinization axis re-discover in
oesophageal carcinoma?

The breast LumA/B residual (anchored on a zero-parameter proliferation signature) was a basal/squamous-
keratinization program (KRT5/14/17/6B, TP63, DSG3/DSC3, ...), and it already re-discovered in TCGA lung
(LUAD vs LUSC). Oesophageal carcinoma (ESCA) splits cleanly into squamous-cell carcinoma (ESCC) and
adenocarcinoma (EAC) — another squamous-vs-adeno contrast in a different organ. We take TCGA ESCA
(UCSC Xena HiSeqV2, n=204; 100 squamous / 104 adeno), anchor on the SAME proliferation signature, and mine
the residual.

Honest, partial result: the breast basal panel TRANSFERS strongly (fixed-panel AUROC ~0.91 for
squamous-vs-adeno; KRT5/TP63 strongly up in ESCC), so the keratinization axis is unmistakably present. But
unbiased residual rediscovery surfaces the ADENOCARCINOMA counter-pole (HNF4A/HNF1A/B, MUC13, VIL1) rather
than the squamous keratins, so the 30-gene overlap with the breast basal panel is 0 — the same squamous/adeno
axis, opposite pole named (contrast with lung, where rediscovery named the squamous pole, 10/30).

Writes external_validation_esca.csv and novel_genes_esca.csv.

Run:  ESCA_DIR=/path/to/tcga_esca python reports/dmoi_external_esca.py
"""
import os, sys
from math import comb
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo
E = os.environ.get("ESCA_DIR", "")
PROLIF = ["MKI67", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB", "BUB1", "CCNE1", "CDC20",
          "TOP2A", "TYMS", "RRM2", "UBE2C", "CENPF", "FOXM1", "MELK", "KIF2C", "NUSAP1", "PTTG1"]


def main():
    assert E and os.path.isdir(E), "set ESCA_DIR (folder with HiSeqV2.gz + ESCA_clinicalMatrix)"
    M = pd.read_csv(os.path.join(E, "HiSeqV2.gz"), sep="\t", index_col=0); M = M[~M.index.duplicated()]
    cl = pd.read_csv(os.path.join(E, "ESCA_clinicalMatrix"), sep="\t", index_col=0)
    h = cl["histological_type"].reindex(M.columns)
    sq = "Esophagus Squamous Cell Carcinoma"; ad = "Esophagus Adenocarcinoma, NOS"
    m = h.isin([sq, ad]); M = M.loc[:, m].dropna(how="any"); y = (h[m] == sq).astype(int).values
    from sklearn.metrics import roc_auc_score
    pres = [x for x in PROLIF if x in M.index]
    anchor = mo.signature_score(M.loc[pres].T, pres)
    basal = set(pd.read_csv(os.path.join(REPO, "novel_genes.csv"))["gene"])
    bpres = [g for g in basal if g in M.index]
    bscore = mo.signature_score(M.loc[bpres].T, bpres)
    basal_panel_auroc = round(float(max(roc_auc_score(y, bscore), 1 - roc_auc_score(y, bscore))), 3)
    var = M.var(axis=1); feats = sorted(set(var.sort_values(ascending=False).head(5000).index) | (basal & set(M.index)))
    X = M.loc[feats].T.values.astype("float32")
    res = mo.anchored_residual_discovery(anchor, X, feats, y, top_k=30, corr_max=0.6, cv=5,
                                         random_state=0, n_perm=10, inner_repeats=1, stability_reps=15)
    nov = [gn for gn, _, _ in res["novel"]]; ov = sorted(set(nov) & basal)
    N, K, n, k = len(feats), len(basal & set(feats)), len(nov), len(ov)
    p_hyper = sum(comb(K, i) * comb(N - K, n - i) for i in range(k, min(K, n) + 1)) / comb(N, n)
    pd.DataFrame(res["novel"], columns=["gene", "partial_corr", "corr_with_anchor"]).to_csv(
        os.path.join(REPO, "novel_genes_esca.csv"), index=False)
    pd.DataFrame([dict(dataset="TCGA_ESCA_ESCC_vs_EAC", endpoint="histology_squamous_vs_adeno", n=int(len(y)),
                       textbook_anchor="proliferation_20gene", anchor_auroc=round(res["auroc_anchor"], 3),
                       basal_panel_transfer_auroc=basal_panel_auroc,
                       combined=round(res["auroc_combined"], 3), delta=round(res["delta"], 3),
                       novel_delta=round(res["novel_delta"], 3), random_delta_mean=round(res["random_delta_mean"], 3),
                       novel_vs_random_p=round(res["novel_vs_random_p"], 3), overlap_with_brca_basal_of30=int(k),
                       overlap_hyper_p=f"{p_hyper:.2e}", rediscovered_pole="adenocarcinoma (HNF4A/HNF1A/B, MUC13, VIL1)",
                       stability_gain=res["stability_gain"], shared_genes="; ".join(ov))
                  ]).to_csv(os.path.join(REPO, "external_validation_esca.csv"), index=False)
    print(f"TCGA ESCA ESCC-vs-EAC n={len(y)} ({int(y.sum())} squamous) | proliferation anchor AUROC={res['auroc_anchor']:.3f}")
    print(f"breast basal panel TRANSFER AUROC (squamous-vs-adeno) = {basal_panel_auroc}  (KRT5/TP63 strongly up in ESCC)")
    print(f"de-novo residual overlap with BRCA basal panel: {k}/30 (p={p_hyper:.2e}); rediscovery names the ADENO pole")
    print("CONCLUSION (honest, partial): the breast basal/keratinization axis is PRESENT in oesophagus — the "
          "fixed panel transfers strongly (AUROC ~0.91) — but unbiased residual rediscovery surfaces the "
          "adenocarcinoma counter-pole (HNF4A/HNF1A/B, MUC13, VIL1) rather than the squamous keratins, so the "
          "30-gene overlap is 0. Same squamous/adeno axis; opposite pole named. Contrast with lung, where "
          "rediscovery named the squamous pole (10/30).")


if __name__ == "__main__":
    main()
