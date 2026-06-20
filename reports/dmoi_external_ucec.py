#!/usr/bin/env python3
"""Cross-cancer validation #7: breast basal/keratinization axis in endometrial cancer (UCEC).

BACKGROUND — SECOND NEGATIVE CONTROL
--------------------------------------
Cross-cancer series:
  #1  Lung  LUAD/LUSC     : overlap 10/30, p=3e-16  → squamous pole        (AUROC ~0.96)
  #2  HNSC  tissue-indep  : panel AUROC 0.962                               (AUROC  0.962)
  #3  ESCA  ESCC vs EAC   : panel AUROC 0.913, 0/30  → adeno counter-pole
  #4  BLCA  Basal/Luminal : panel AUROC 0.967, 1/30  → luminal counter-pole
  #5  CESC  Sq vs Adeno   : panel AUROC 0.938, 6/30  → squamous pole
  #6  STAD  Int vs Diff   : panel AUROC 0.517 [NEG]  → immune infiltration  (no squamous)
  #7  UCEC  (this script) : SECOND NEGATIVE CONTROL                         (no squamous)

UCEC (uterine corpus endometrial carcinoma) — BOTH poles are adenocarcinoma:
  - Endometrioid endometrial adenocarcinoma (y=0): ER+, low-grade, microsatellite stable/MSI
  - Serous endometrial adenocarcinoma        (y=1): TP53-mutant, HER2-amplified, aggressive

PREDICTION:
  - Breast basal panel TRANSFER AUROC ≈ 0.5 (no squamous component → should not transfer)
  - Residual overlap ≈ 0/30
  - Residual will name endometrioid-specific (FOXA2, CDX2, PTEN-pathway) or
    serous-specific markers (TP53-target genes, CCNE1, MUC16/CA-125)
  - Second data point confirming adeno-vs-adeno → axis specificity

TCGA UCEC (UCSC Xena HiSeqV2, n=562):
  y=1  Serous      : "Serous endometrial adenocarcinoma"                 n=121
  y=0  Endometrioid: "Endometrioid endometrial adenocarcinoma"           n=441
  Excluded         : "Mixed serous and endometrioid"                     n= 22

RUN
---
  python reports/dmoi_external_ucec.py
  UCEC_DIR=/path/to/cache python reports/dmoi_external_ucec.py
"""

import os, sys, urllib.request
from math import comb
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo

XENA_EXPR = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
    "TCGA.UCEC.sampleMap/HiSeqV2.gz"
)
XENA_CLIN = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
    "TCGA.UCEC.sampleMap/UCEC_clinicalMatrix"
)

SER_CLASS  = "Serous endometrial adenocarcinoma"           # y=1
ENDO_CLASS = "Endometrioid endometrial adenocarcinoma"     # y=0
# "Mixed serous and endometrioid" excluded

PROLIF = [
    "MKI67", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB",
    "BUB1", "CCNE1", "CDC20", "TOP2A", "TYMS", "RRM2", "UBE2C",
    "CENPF", "FOXM1", "MELK", "KIF2C", "NUSAP1", "PTTG1",
]


def _fetch(url, dest):
    if not os.path.exists(dest):
        print(f"  downloading {os.path.basename(dest)} …", flush=True)
        urllib.request.urlretrieve(url, dest)


def load_data():
    cache = os.environ.get("UCEC_DIR", os.getcwd())
    expr_path = os.path.join(cache, "UCEC_HiSeqV2.gz")
    clin_path  = os.path.join(cache, "UCEC_clinicalMatrix")
    _fetch(XENA_EXPR, expr_path)
    _fetch(XENA_CLIN, clin_path)

    M  = pd.read_csv(expr_path, sep="\t", index_col=0)
    M  = M[~M.index.duplicated()]
    cl = pd.read_csv(clin_path, sep="\t", index_col=0)

    h    = cl["histological_type"].reindex(M.columns)
    keep = h.isin([SER_CLASS, ENDO_CLASS])
    M    = M.loc[:, keep].dropna(how="any")
    y    = (h[keep] == SER_CLASS).astype(int).values

    n_ser = int(y.sum()); n_end = int((1-y).sum())
    print(f"  serous={n_ser}  endometrioid={n_end}  total={len(y)}", flush=True)
    return M, y


def main():
    print("== UCEC cross-cancer validation #7 — SECOND NEGATIVE CONTROL ==", flush=True)
    print("   [NO squamous: serous vs endometrioid, both adenocarcinoma]", flush=True)
    M, y = load_data()

    pres   = [x for x in PROLIF if x in M.index]
    anchor = mo.signature_score(M.loc[pres].T, pres)

    basal  = set(pd.read_csv(os.path.join(REPO, "novel_genes.csv"))["gene"])
    bpres  = [g for g in basal if g in M.index]
    bscore = mo.signature_score(M.loc[bpres].T, bpres)
    raw    = roc_auc_score(y, bscore)
    basal_panel_auroc = round(float(max(raw, 1 - raw)), 3)

    var   = M.var(axis=1)
    feats = sorted(
        set(var.sort_values(ascending=False).head(5000).index) | (basal & set(M.index))
    )
    X = M.loc[feats].T.values.astype("float32")
    res = mo.anchored_residual_discovery(
        anchor, X, feats, y,
        top_k=30, corr_max=0.6, cv=5,
        random_state=0, n_perm=10, inner_repeats=1, stability_reps=15,
    )

    nov  = [gn for gn, _, _ in res["novel"]]
    ov   = sorted(set(nov) & basal)
    N, K, n_nov, k = len(feats), len(basal & set(feats)), len(nov), len(ov)
    p_hyper = (
        sum(comb(K, i) * comb(N - K, n_nov - i) for i in range(k, min(K, n_nov) + 1))
        / comb(N, n_nov)
    )

    pd.DataFrame(res["novel"], columns=["gene", "partial_corr", "corr_with_anchor"]).to_csv(
        os.path.join(REPO, "novel_genes_ucec.csv"), index=False
    )
    pd.DataFrame([dict(
        dataset                      = "TCGA_UCEC_Serous_vs_Endometrioid",
        endpoint                     = "histology_serous_vs_endometrioid",
        n                            = int(len(y)),
        n_serous                     = int(y.sum()),
        n_endometrioid               = int((1-y).sum()),
        textbook_anchor              = "proliferation_20gene",
        anchor_auroc                 = round(res["auroc_anchor"], 3),
        basal_panel_transfer_auroc   = basal_panel_auroc,
        combined                     = round(res["auroc_combined"], 3),
        delta                        = round(res["delta"], 3),
        novel_delta                  = round(res["novel_delta"], 3),
        random_delta_mean            = round(res["random_delta_mean"], 3),
        novel_vs_random_p            = round(res["novel_vs_random_p"], 3),
        overlap_with_brca_basal_of30 = int(k),
        overlap_hyper_p              = f"{p_hyper:.2e}",
        rediscovered_pole            = "?",
        stability_gain               = res["stability_gain"],
        shared_genes                 = "; ".join(ov),
        note                         = "NEGATIVE CONTROL — no squamous component",
    )]).to_csv(os.path.join(REPO, "external_validation_ucec.csv"), index=False)

    print(f"\n  anchor AUROC (prolif vs serous/endometrioid): {res['auroc_anchor']:.3f}")
    print(f"  breast basal panel TRANSFER AUROC:            {basal_panel_auroc}")
    print(f"  residual overlap with BRCA basal panel:       {k}/30  (p={p_hyper:.2e})")
    if ov: print(f"  shared genes: {', '.join(ov)}")
    print(f"  top 15 residual genes: {', '.join(nov[:15])}")

    print()
    if basal_panel_auroc < 0.65:
        verdict = f"NEGATIVE (AUROC={basal_panel_auroc} ≈ chance) — confirms adeno-vs-adeno specificity"
    elif basal_panel_auroc < 0.80:
        verdict = f"WEAK ({basal_panel_auroc}) — investigate squamous elements"
    else:
        verdict = f"UNEXPECTED HIGH ({basal_panel_auroc}) — INVESTIGATE"
    print(f"  *** VERDICT: {verdict}")

    print("\n" + "=" * 62)
    print("CROSS-CANCER SERIES SUMMARY")
    print("=" * 62)
    print("  Squamous-containing:")
    print("    #1 Lung  LUAD/LUSC   : overlap 10/30, p=3.2e-16  → squamous pole")
    print("    #2 HNSC  tissue-indep: panel AUROC 0.962")
    print("    #3 ESCA  ESCC/EAC    : panel AUROC 0.913, 0/30   → adeno pole")
    print("    #4 BLCA  Basal/Lum   : panel AUROC 0.967, 1/30   → luminal pole")
    print("    #5 CESC  Sq/Adeno    : panel AUROC 0.938, 6/30   → squamous pole")
    print("  Negative controls (no squamous):")
    print("    #6 STAD  Int/Diff    : panel AUROC 0.517 [NEG]   → immune infiltration")
    print(f"    #7 UCEC  Ser/Endo    : panel AUROC {basal_panel_auroc}   "
          f"overlap {k}/30 [NEG CTRL]")


if __name__ == "__main__":
    main()
