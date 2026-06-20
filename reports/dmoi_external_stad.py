#!/usr/bin/env python3
"""Cross-cancer validation #6: breast basal/keratinization axis in gastric cancer (STAD).

BACKGROUND — NEGATIVE CONTROL EXPERIMENT
-----------------------------------------
Cross-cancer series:
  #1  Lung  LUAD/LUSC     : overlap 10/30, p=3e-16  → squamous pole (HIGH transfer)
  #2  HNSC  tissue-indep  : panel AUROC 0.962        (HIGH transfer)
  #3  ESCA  ESCC vs EAC   : panel AUROC 0.913, 0/30  → adeno counter-pole
  #4  BLCA  Basal/Luminal : panel AUROC 0.967, 1/30  → luminal/urothelial counter-pole
  #5  CESC  Sq vs Adeno   : panel AUROC 0.938, 6/30  → squamous pole (HPV-amplified)
  #6  STAD  (this script) : NEGATIVE CONTROL — no squamous component

STAD (stomach adenocarcinoma) is ENTIRELY adenocarcinoma.
It has TWO biologically distinct subtypes (Lauren classification analogue):
  - Intestinal-type (y=1): tubular/mucinous/papillary; CDX2+, MUC2+, HER2-amp
  - Diffuse-type    (y=0): signet ring cells; CDH1 loss/mutation, RHOA mutation

PREDICTION (squamous-axis specificity hypothesis):
  - Breast basal panel TRANSFER AUROC ≈ 0.5 (panel should NOT transfer — no squamous axis)
  - Residual overlap ≈ 0/30
  - Residual will name intestinal-type markers (CDX2, MUC2, VIL1) or diffuse markers (CDH1)
  - If panel AUROC is truly ≈ 0.5 → PROVES SPECIFICITY of the pan-epithelial keratinization
    axis to squamous-containing comparisons (not generic epithelial differentiation)

HISTOLOGICAL TYPES USED (from TCGA STAD clinical matrix):
  y=1 Intestinal: "Stomach, Intestinal Adenocarcinoma, Tubular Type"    n=85
                + "Stomach, Intestinal Adenocarcinoma, Mucinous Type"   n=23
                + "Stomach, Intestinal Adenocarcinoma, Papillary Type"  n=9
                TOTAL: n=117

  y=0 Diffuse:   "Stomach, Adenocarcinoma, Diffuse Type"                n=80
                + "Stomach Adenocarcinoma, Signet Ring Type"             n=13
                TOTAL: n=93

  EXCLUDED:      "NOS" types (n=208+90=298, too ambiguous)

RUN
---
  python reports/dmoi_external_stad.py
  STAD_DIR=/path/to/cache python reports/dmoi_external_stad.py
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
    "TCGA.STAD.sampleMap/HiSeqV2.gz"
)
XENA_CLIN = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
    "TCGA.STAD.sampleMap/STAD_clinicalMatrix"
)

# ── Class definitions ──────────────────────────────────────────────────────────
INT_CLASSES = {                        # y=1  intestinal-type
    "Stomach, Intestinal Adenocarcinoma, Tubular Type",
    "Stomach, Intestinal Adenocarcinoma, Mucinous Type",
    "Stomach, Intestinal Adenocarcinoma, Papillary Type",
}
DIFF_CLASSES = {                       # y=0  diffuse-type
    "Stomach, Adenocarcinoma, Diffuse Type",
    "Stomach Adenocarcinoma, Signet Ring Type",
}
# "NOS" types excluded — too ambiguous

# ── Anchor (constant across series) ───────────────────────────────────────────
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
    cache = os.environ.get("STAD_DIR", os.getcwd())
    expr_path = os.path.join(cache, "STAD_HiSeqV2.gz")
    clin_path  = os.path.join(cache, "STAD_clinicalMatrix")
    _fetch(XENA_EXPR, expr_path)
    _fetch(XENA_CLIN, clin_path)

    M  = pd.read_csv(expr_path, sep="\t", index_col=0)
    M  = M[~M.index.duplicated()]
    cl = pd.read_csv(clin_path, sep="\t", index_col=0)

    h    = cl["histological_type"].reindex(M.columns)
    keep = h.isin(INT_CLASSES | DIFF_CLASSES)
    M    = M.loc[:, keep].dropna(how="any")
    y    = (h[keep].isin(INT_CLASSES)).astype(int).values

    n_int = int(y.sum()); n_diff = int((1-y).sum())
    print(f"  intestinal={n_int}  diffuse/signet={n_diff}  total={len(y)}", flush=True)
    return M, y


def main():
    print("== STAD cross-cancer validation #6 — NEGATIVE CONTROL ==", flush=True)
    print("   [NO squamous component: testing pan-epithelial axis SPECIFICITY]", flush=True)
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
        os.path.join(REPO, "novel_genes_stad.csv"), index=False
    )

    # infer pole label
    pole = "intestinal adenocarcinoma (CDX2/MUC2/VIL1 expected)"
    if any(g in nov[:10] for g in ("CDH1", "RHOA", "RHO")):
        pole = "diffuse/signet-ring counter-pole"

    pd.DataFrame([dict(
        dataset                      = "TCGA_STAD_Intestinal_vs_Diffuse",
        endpoint                     = "histology_intestinal_vs_diffuse",
        n                            = int(len(y)),
        n_intestinal                 = int(y.sum()),
        n_diffuse                    = int((1-y).sum()),
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
        rediscovered_pole            = pole,
        stability_gain               = res["stability_gain"],
        shared_genes                 = "; ".join(ov),
        note                         = "NEGATIVE CONTROL — no squamous component",
    )]).to_csv(os.path.join(REPO, "external_validation_stad.csv"), index=False)

    # print
    print(f"\n  anchor AUROC (prolif vs intestinal/diffuse): {res['auroc_anchor']:.3f}")
    print(f"  breast basal panel TRANSFER AUROC:           {basal_panel_auroc}")
    print(f"  residual overlap with BRCA basal panel:      {k}/30  (p={p_hyper:.2e})")
    if ov:
        print(f"  shared genes: {', '.join(ov)}")
    print(f"  top 15 residual genes: {', '.join(nov[:15])}")

    # verdict
    print()
    if basal_panel_auroc < 0.65:
        print(f"  *** VERDICT: NEGATIVE (panel AUROC={basal_panel_auroc} ≈ chance)")
        print(f"      → Confirms specificity of pan-epithelial axis to squamous biology")
        print(f"      → Panel DOES NOT transfer to adeno-vs-adeno comparisons")
    elif basal_panel_auroc < 0.80:
        print(f"  *** VERDICT: WEAK PARTIAL (panel AUROC={basal_panel_auroc})")
        print(f"      → Unexpected: some keratinization signal in intestinal vs diffuse")
        print(f"      → Investigate: intestinal-type may have squamous metaplasia elements")
    else:
        print(f"  *** UNEXPECTED: panel AUROC={basal_panel_auroc} — INVESTIGATE")
        print(f"      → If high overlap → re-examine endpoint purity")

    # series summary
    print("\n" + "=" * 62)
    print("CROSS-CANCER SERIES SUMMARY")
    print("=" * 62)
    print("  #1 Lung  LUAD/LUSC     : overlap 10/30, p=3.2e-16 → squamous pole")
    print("  #2 HNSC  tissue-indep  : panel AUROC 0.962")
    print("  #3 ESCA  ESCC vs EAC   : panel AUROC 0.913, 0/30 → adeno pole")
    print("  #4 BLCA  Basal/Luminal : panel AUROC 0.967, 1/30 → luminal pole")
    print("  #5 CESC  Sq vs Adeno   : panel AUROC 0.938, 6/30 → squamous pole")
    print(f"  #6 STAD  Int vs Diff   : panel AUROC {basal_panel_auroc}  "
          f"overlap {k}/30 (p={p_hyper:.2e}) ← NEGATIVE CONTROL")


if __name__ == "__main__":
    main()
