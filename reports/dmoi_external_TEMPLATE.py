#!/usr/bin/env python3
"""Cross-cancer validation TEMPLATE — copy and fill in the ### CONFIG ### section.

USAGE
-----
1. Copy this file:  cp reports/dmoi_external_TEMPLATE.py reports/dmoi_external_<cancer>.py
2. Fill every line marked  ### CONFIG ###
3. Run:  <CANCER>_DIR=/path/to/xena_folder python reports/dmoi_external_<cancer>.py
         (or leave <CANCER>_DIR empty to trigger Xena auto-download, if _get() is implemented)

WHAT THIS SCRIPT DOES
---------------------
Anchored residual discovery cross-cancer validation:

  Step 1.  Load expression matrix M (genes × samples) and binary label y.
  Step 2.  Anchor score = mean of 20-gene proliferation signature per sample.
  Step 3.  PANEL TRANSFER: score each sample with the breast basal panel (novel_genes.csv);
           compute AUROC for distinguishing y with that fixed panel.
  Step 4.  RESIDUAL DISCOVERY: anchored_residual_discovery() — partial-correlates out the
           anchor, then mines residual features that improve y-prediction.
  Step 5.  Overlap with breast basal panel (hypergeometric p).
  Step 6.  Write <CANCER>_validation.csv and novel_genes_<cancer>.csv.

Cross-cancer series so far:
  #1  LUAD vs LUSC  (lung):       overlap 10/30, p=3e-16   → squamous pole re-discovered
  #2  HNSC tissue-independence:   AUROC 0.962 squamous panel
  #3  ESCA (ESCC vs EAC):         overlap  0/30, p=1.0     → adeno pole surfaced
      → same squamous/adeno axis; which pole is named differs by starting tissue
"""

import os, sys, urllib.request
from math import comb
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo
from sklearn.metrics import roc_auc_score

# ─────────────────────────────────────────────────────────────────────
# ### CONFIG ###  — fill every value in this section
# ─────────────────────────────────────────────────────────────────────

CANCER_TAG     = "BLCA"               # ### CONFIG ### short tag used in filenames
DATASET_LABEL  = "TCGA_BLCA_??_vs_??"# ### CONFIG ### value for CSV 'dataset' column
ENDPOINT_LABEL = "??"                 # ### CONFIG ### e.g. "subtype_luminal_vs_basal"

# Xena URL template — replace {CANCER_TAG} if the Xena cohort name differs
XENA_EXPR_URL  = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
    f"TCGA.{CANCER_TAG}.sampleMap/HiSeqV2.gz"
)
XENA_CLIN_URL  = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
    f"TCGA.{CANCER_TAG}.sampleMap/{CANCER_TAG}_clinicalMatrix"
)

# Environment variable the user sets to a local folder (bypasses download)
ENV_VAR        = f"{CANCER_TAG}_DIR"  # ### CONFIG ### e.g. "BLCA_DIR"

# Clinical matrix column and the two values that define positive / negative class
CLINICAL_COL   = "histological_type"  # ### CONFIG ### column in clinicalMatrix
POS_CLASS      = "??"                 # ### CONFIG ### y=1 (e.g. "Transitional Cell Carcinoma")
NEG_CLASS      = "??"                 # ### CONFIG ### y=0

# Human-readable description of what pole you EXPECT residual discovery to name.
# Fill after running — leave as "?" for the first run.
EXPECTED_POLE  = "?"                  # ### CONFIG ###

# ─────────────────────────────────────────────────────────────────────
# Anchor gene panel (DO NOT CHANGE — same across all cross-cancer tests)
# ─────────────────────────────────────────────────────────────────────
PROLIF = [
    "MKI67", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB",
    "BUB1", "CCNE1", "CDC20", "TOP2A", "TYMS", "RRM2", "UBE2C",
    "CENPF", "FOXM1", "MELK", "KIF2C", "NUSAP1", "PTTG1",
]

OUTPUT_CSV       = os.path.join(REPO, f"external_validation_{CANCER_TAG.lower()}.csv")
NOVEL_GENES_CSV  = os.path.join(REPO, f"novel_genes_{CANCER_TAG.lower()}.csv")


# ─────────────────────────────────────────────────────────────────────
# Data loading helpers
# ─────────────────────────────────────────────────────────────────────

def _cache_dir():
    d = os.environ.get(ENV_VAR, "")
    if d and os.path.isdir(d):
        return d
    # auto-download into cwd if ENV_VAR not set
    return os.getcwd()


def _fetch(url, dest):
    if not os.path.exists(dest):
        print(f"  downloading {os.path.basename(dest)} …")
        urllib.request.urlretrieve(url, dest)


def load_data():
    """Return (M: pd.DataFrame genes×samples, y: np.ndarray int) with y=1 for POS_CLASS."""
    d = _cache_dir()
    expr_path = os.path.join(d, f"{CANCER_TAG}_HiSeqV2.gz")
    clin_path  = os.path.join(d, f"{CANCER_TAG}_clinicalMatrix")

    _fetch(XENA_EXPR_URL, expr_path)
    _fetch(XENA_CLIN_URL, clin_path)

    M  = pd.read_csv(expr_path, sep="\t", index_col=0)
    M  = M[~M.index.duplicated()]
    cl = pd.read_csv(clin_path, sep="\t", index_col=0)

    labels = cl[CLINICAL_COL].reindex(M.columns)
    keep   = labels.isin([POS_CLASS, NEG_CLASS])
    M      = M.loc[:, keep].dropna(how="any")
    y      = (labels[keep] == POS_CLASS).astype(int).values
    return M, y


# ─────────────────────────────────────────────────────────────────────
# Main analysis
# ─────────────────────────────────────────────────────────────────────

def main():
    print(f"== {CANCER_TAG} cross-cancer validation ==")
    print(f"   POS={POS_CLASS!r}  NEG={NEG_CLASS!r}")

    M, y = load_data()
    n_pos, n_neg = int(y.sum()), int((1 - y).sum())
    print(f"   n={len(y)}  ({n_pos} pos / {n_neg} neg)")

    # --- anchor ---
    pres   = [x for x in PROLIF if x in M.index]
    anchor = mo.signature_score(M.loc[pres].T, pres)

    # --- panel transfer (breast basal panel → this cancer) ---
    basal  = set(pd.read_csv(os.path.join(REPO, "novel_genes.csv"))["gene"])
    bpres  = [g for g in basal if g in M.index]
    bscore = mo.signature_score(M.loc[bpres].T, bpres)
    raw_auroc = roc_auc_score(y, bscore)
    basal_panel_auroc = round(float(max(raw_auroc, 1 - raw_auroc)), 3)

    # --- residual discovery ---
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

    # --- overlap statistics ---
    nov  = [gn for gn, _, _ in res["novel"]]
    ov   = sorted(set(nov) & basal)
    N, K, n_nov, k = len(feats), len(basal & set(feats)), len(nov), len(ov)
    p_hyper = (
        sum(comb(K, i) * comb(N - K, n_nov - i) for i in range(k, min(K, n_nov) + 1))
        / comb(N, n_nov)
    )

    # --- save ---
    pd.DataFrame(res["novel"], columns=["gene", "partial_corr", "corr_with_anchor"]).to_csv(
        NOVEL_GENES_CSV, index=False
    )
    pd.DataFrame([dict(
        dataset                  = DATASET_LABEL,
        endpoint                 = ENDPOINT_LABEL,
        n                        = int(len(y)),
        n_pos                    = n_pos,
        n_neg                    = n_neg,
        textbook_anchor          = "proliferation_20gene",
        anchor_auroc             = round(res["auroc_anchor"], 3),
        basal_panel_transfer_auroc = basal_panel_auroc,
        combined                 = round(res["auroc_combined"], 3),
        delta                    = round(res["delta"], 3),
        novel_delta              = round(res["novel_delta"], 3),
        random_delta_mean        = round(res["random_delta_mean"], 3),
        novel_vs_random_p        = round(res["novel_vs_random_p"], 3),
        overlap_with_brca_basal_of30 = int(k),
        overlap_hyper_p          = f"{p_hyper:.2e}",
        rediscovered_pole        = EXPECTED_POLE,
        stability_gain           = res["stability_gain"],
        shared_genes             = "; ".join(ov),
    )]).to_csv(OUTPUT_CSV, index=False)

    # --- print ---
    print(f"\n  anchor AUROC (prolif vs {ENDPOINT_LABEL}): {res['auroc_anchor']:.3f}")
    print(f"  breast basal panel TRANSFER AUROC:          {basal_panel_auroc}")
    print(f"  residual overlap with BRCA basal panel:     {k}/30  (p={p_hyper:.2e})")
    if ov:
        print(f"  shared genes: {', '.join(ov)}")
    print(f"\n  top residual genes: {', '.join(nov[:10])}")
    print(f"\n  → fill EXPECTED_POLE above based on what you see in top residual genes")
    print(f"  saved: {OUTPUT_CSV}")
    print(f"         {NOVEL_GENES_CSV}")


if __name__ == "__main__":
    main()
