#!/usr/bin/env python3
"""Cross-cancer validation #5: breast basal/keratinization axis in cervical carcinoma (CESC).

BACKGROUND
----------
Cross-cancer series:
  #1  Lung LUAD/LUSC     : overlap 10/30, p=3e-16   → squamous pole (same pole)
  #2  HNSC tissue-indep  : panel AUROC 0.962
  #3  ESCA ESCC vs EAC   : panel AUROC 0.913, overlap  0/30 → adeno pole (counter-pole)
  #4  BLCA Basal/Luminal : panel AUROC 0.967, overlap  1/30 → luminal/urothelial (counter-pole)
  #5  CESC (this script) : ?

CESC splits into:
  - Cervical Squamous Cell Carcinoma (CSCC): HPV+ squamous, TP63/KRT5/KRT14 high
    → PREDICTION: breast basal panel scores HIGH
  - Endocervical Adenocarcinoma (EAC): glandular, PAX8/HNF1A/MUC16 high
    → PREDICTION: residual will name adeno counter-pole (as in ESCA) OR squamous pole (as in lung)

Key question: does cervical squamous carcinoma share the keratinization axis, and what does
de-novo residual discovery surface — the squamous pole (→ lung pattern) or the adeno counter-pole
(→ ESCA/BLCA pattern)?

TCGA CESC (UCSC Xena HiSeqV2, n ≈ 305):
  - Squamous (y=1): "Cervical Squamous Cell Carcinoma"           n=257
  - Adeno    (y=0): four endocervical adenocarcinoma subtypes    n= 48
  - Excluded      : "Adenosquamous"                              n=  7

RUN
---
  # auto-download from Xena:
  python reports/dmoi_external_cesc.py

  # use local cache:
  CESC_DIR=/path/to/cache python reports/dmoi_external_cesc.py
"""

import os, sys, urllib.request
from math import comb
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo

# ── Xena URLs ──────────────────────────────────────────────────────────────────
XENA_EXPR = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
    "TCGA.CESC.sampleMap/HiSeqV2.gz"
)
XENA_CLIN = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
    "TCGA.CESC.sampleMap/CESC_clinicalMatrix"
)

# ── Histological class labels ───────────────────────────────────────────────────
SQ_CLASS = "Cervical Squamous Cell Carcinoma"          # y=1
AD_CLASSES = {                                          # y=0 (all adeno subtypes)
    "Endocervical Type of Adenocarcinoma",
    "Mucinous Adenocarcinoma of Endocervical Type",
    "Endocervical Adenocarcinoma of the Usual Type",
    "Endometrioid Adenocarcinoma of Endocervix",
}
# "Adenosquamous" excluded

# ── Anchor gene panel (same across ALL cross-cancer tests) ─────────────────────
PROLIF = [
    "MKI67", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB",
    "BUB1", "CCNE1", "CDC20", "TOP2A", "TYMS", "RRM2", "UBE2C",
    "CENPF", "FOXM1", "MELK", "KIF2C", "NUSAP1", "PTTG1",
]


# ── Helpers ────────────────────────────────────────────────────────────────────

def _fetch(url, dest):
    if not os.path.exists(dest):
        print(f"  downloading {os.path.basename(dest)} …", flush=True)
        urllib.request.urlretrieve(url, dest)


def load_data():
    cache = os.environ.get("CESC_DIR", os.getcwd())
    expr_path = os.path.join(cache, "CESC_HiSeqV2.gz")
    clin_path  = os.path.join(cache, "CESC_clinicalMatrix")
    _fetch(XENA_EXPR, expr_path)
    _fetch(XENA_CLIN, clin_path)

    M  = pd.read_csv(expr_path, sep="\t", index_col=0)
    M  = M[~M.index.duplicated()]
    cl = pd.read_csv(clin_path, sep="\t", index_col=0)

    h    = cl["histological_type"].reindex(M.columns)
    keep = h.isin({SQ_CLASS} | AD_CLASSES)
    M    = M.loc[:, keep].dropna(how="any")
    y    = (h[keep] == SQ_CLASS).astype(int).values

    n_sq = int(y.sum()); n_ad = int((1-y).sum())
    print(f"  squamous={n_sq}  adeno={n_ad}  total={len(y)}", flush=True)
    return M, y


# ── Core analysis ──────────────────────────────────────────────────────────────

def main():
    print("== CESC cross-cancer validation #5 ==", flush=True)
    M, y = load_data()

    # anchor
    pres   = [x for x in PROLIF if x in M.index]
    anchor = mo.signature_score(M.loc[pres].T, pres)

    # panel transfer
    basal  = set(pd.read_csv(os.path.join(REPO, "novel_genes.csv"))["gene"])
    bpres  = [g for g in basal if g in M.index]
    bscore = mo.signature_score(M.loc[bpres].T, bpres)
    raw    = roc_auc_score(y, bscore)
    basal_panel_auroc = round(float(max(raw, 1 - raw)), 3)

    # residual discovery
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

    # overlap
    nov  = [gn for gn, _, _ in res["novel"]]
    ov   = sorted(set(nov) & basal)
    N, K, n_nov, k = len(feats), len(basal & set(feats)), len(nov), len(ov)
    p_hyper = (
        sum(comb(K, i) * comb(N - K, n_nov - i) for i in range(k, min(K, n_nov) + 1))
        / comb(N, n_nov)
    )

    # save
    pd.DataFrame(res["novel"], columns=["gene", "partial_corr", "corr_with_anchor"]).to_csv(
        os.path.join(REPO, "novel_genes_cesc.csv"), index=False
    )
    pd.DataFrame([dict(
        dataset                      = "TCGA_CESC_Squamous_vs_Adeno",
        endpoint                     = "histology_squamous_vs_adeno",
        n                            = int(len(y)),
        n_squamous                   = int(y.sum()),
        n_adeno                      = int((1-y).sum()),
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
        rediscovered_pole            = "?",   # fill after inspecting nov
        stability_gain               = res["stability_gain"],
        shared_genes                 = "; ".join(ov),
    )]).to_csv(os.path.join(REPO, "external_validation_cesc.csv"), index=False)

    # print
    print(f"  anchor AUROC (prolif vs squamous/adeno): {res['auroc_anchor']:.3f}")
    print(f"  breast basal panel TRANSFER AUROC:       {basal_panel_auroc}")
    print(f"  residual overlap with BRCA basal panel:  {k}/30  (p={p_hyper:.2e})")
    if ov:
        print(f"  shared genes: {', '.join(ov)}")
    print(f"  top 15 residual genes: {', '.join(nov[:15])}")

    # series summary
    print("\n" + "=" * 62)
    print("CROSS-CANCER SERIES SUMMARY")
    print("=" * 62)
    print("  #1 Lung  LUAD/LUSC     : overlap 10/30, p=3.2e-16  → squamous pole")
    print("  #2 HNSC  tissue-indep  : panel AUROC 0.962")
    print("  #3 ESCA  ESCC vs EAC   : panel AUROC 0.913, overlap  0/30 → adeno pole")
    print("  #4 BLCA  Basal/Luminal : panel AUROC 0.967, overlap  1/30 → luminal/urothelial")
    print(f"  #5 CESC  Sq vs Adeno  : panel AUROC {basal_panel_auroc}, "
          f"overlap {k}/30 (p={p_hyper:.2e})")
    print()
    print("  PATTERN KEY:")
    print("    overlap >> 0 → residual names squamous pole (lung pattern)")
    print("    overlap ≈  0 + panel AUROC high → residual names adeno counter-pole (ESCA/BLCA pattern)")
    print()
    print(f"  → inspect top residual genes to determine pole: {', '.join(nov[:8])}")


if __name__ == "__main__":
    main()
