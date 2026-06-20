#!/usr/bin/env python3
"""Cross-cancer validation #4: breast basal/keratinization axis in bladder carcinoma (BLCA).

BACKGROUND
----------
Cross-cancer series so far:
  #1  LUAD vs LUSC  (lung):       overlap 10/30, p=3e-16   → squamous pole re-discovered
  #2  HNSC tissue-independence:   AUROC 0.962 squamous panel
  #3  ESCA (ESCC vs EAC):         overlap  0/30, p=1.0     → adeno pole surfaced
  #4  BLCA (this script):         ?

Bladder carcinoma has two archetypal molecular subtypes (Robertson et al. 2017 Cell):
  - Basal-squamous   : TP63+, KRT5/14+, squamous differentiation markers
                       → PREDICTION: breast basal panel should TRANSFER (high AUROC)
  - Luminal-papillary: FGFR3+, GATA3+, urothelial differentiation
                       → PREDICTION: residual discovery will name the BASAL pole
                         (KRT5/TP63 etc.), giving FIRST overlap > 0 outside lung

This is the critical test: if BLCA basal-squamous shares the keratinization program
with breast basal AND lung squamous, the axis is PAN-EPITHELIAL, not tissue-specific.

TWO PARTS
─────────
Part A  — Clinical proxy (fully automatic, no extra files needed):
  diagnosis_subtype: Non-Papillary (y=1, ≈ Basal-squamous)
                  vs Papillary     (y=0, ≈ Luminal-papillary)
  n ≈ 431  (293 Non-Papillary / 138 Papillary)
  Writes: external_validation_blca_clinical.csv, novel_genes_blca_clinical.csv

Part B  — Robertson 2017 molecular subtypes (requires annotation TSV):
  Basal-squamous (y=1) vs Luminal-papillary (y=0)
  Set env var BLCA_SUBTYPE_TSV to a TSV with columns [sample, subtype].
  Expected subtype values: 'Basal-squamous', 'Luminal-papillary'
  (Neuronal / Luminal / Luminal-infiltrated rows silently dropped.)
  Writes: external_validation_blca_subtype.csv, novel_genes_blca_subtype.csv

  Download supplementary Table S1 from Robertson et al. 2017:
    https://doi.org/10.1016/j.cell.2017.09.007  (Table S1, sheet "TCGA BLCA Subtypes")

RUN
───
  # Part A only (Xena auto-download):
  python reports/dmoi_external_blca.py

  # Both parts:
  BLCA_SUBTYPE_TSV=/path/to/robertson2017_subtypes.tsv \\
    python reports/dmoi_external_blca.py

  # Use local cache to skip re-download:
  BLCA_DIR=/path/to/cache BLCA_SUBTYPE_TSV=... python reports/dmoi_external_blca.py
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
    "TCGA.BLCA.sampleMap/HiSeqV2.gz"
)
XENA_CLIN = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
    "TCGA.BLCA.sampleMap/BLCA_clinicalMatrix"
)

# ── Anchor gene panel (identical across ALL cross-cancer tests — do not change) ─
PROLIF = [
    "MKI67", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB",
    "BUB1", "CCNE1", "CDC20", "TOP2A", "TYMS", "RRM2", "UBE2C",
    "CENPF", "FOXM1", "MELK", "KIF2C", "NUSAP1", "PTTG1",
]


# ── I/O helpers ────────────────────────────────────────────────────────────────

def _fetch(url, dest):
    if not os.path.exists(dest):
        print(f"  downloading {os.path.basename(dest)} …", flush=True)
        urllib.request.urlretrieve(url, dest)


def _load_xena(cache_dir):
    """Return (M: DataFrame genes×samples, cl: DataFrame clinical)."""
    expr_path = os.path.join(cache_dir, "BLCA_HiSeqV2.gz")
    clin_path  = os.path.join(cache_dir, "BLCA_clinicalMatrix")
    _fetch(XENA_EXPR, expr_path)
    _fetch(XENA_CLIN, clin_path)
    M  = pd.read_csv(expr_path, sep="\t", index_col=0)
    M  = M[~M.index.duplicated()]
    cl = pd.read_csv(clin_path, sep="\t", index_col=0)
    return M, cl


# ── Core analysis (shared by Part A and Part B) ────────────────────────────────

def _run(M, y, dataset_label, endpoint_label, expected_pole, out_csv, novel_csv):
    """Anchored residual discovery + overlap stats + file output."""
    n_pos, n_neg = int(y.sum()), int((1 - y).sum())
    print(f"  n={len(y)}  (pos={n_pos} / neg={n_neg})", flush=True)

    # anchor score
    pres   = [x for x in PROLIF if x in M.index]
    anchor = mo.signature_score(M.loc[pres].T, pres)

    # breast basal panel transfer
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

    # overlap statistics
    nov  = [gn for gn, _, _ in res["novel"]]
    ov   = sorted(set(nov) & basal)
    N, K, n_nov, k = len(feats), len(basal & set(feats)), len(nov), len(ov)
    p_hyper = (
        sum(comb(K, i) * comb(N - K, n_nov - i) for i in range(k, min(K, n_nov) + 1))
        / comb(N, n_nov)
    )

    # save
    pd.DataFrame(res["novel"], columns=["gene", "partial_corr", "corr_with_anchor"]).to_csv(
        novel_csv, index=False
    )
    pd.DataFrame([dict(
        dataset                      = dataset_label,
        endpoint                     = endpoint_label,
        n                            = int(len(y)),
        n_pos                        = n_pos,
        n_neg                        = n_neg,
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
        rediscovered_pole            = expected_pole,
        stability_gain               = res["stability_gain"],
        shared_genes                 = "; ".join(ov),
    )]).to_csv(out_csv, index=False)

    # print
    print(f"  anchor AUROC (prolif vs {endpoint_label}): {res['auroc_anchor']:.3f}")
    print(f"  breast basal panel TRANSFER AUROC:          {basal_panel_auroc}")
    print(f"  residual overlap with BRCA basal panel:     {k}/30  (p={p_hyper:.2e})")
    if ov:
        print(f"  shared genes: {', '.join(ov)}")
    print(f"  top 10 residual genes: {', '.join(nov[:10])}")
    print(f"  → update 'expected_pole' in out_csv once you inspect residual genes")
    print(f"  saved: {out_csv}")

    return dict(basal_panel_auroc=basal_panel_auroc, overlap=k, p_hyper=p_hyper, novel=nov)


# ── Part A: clinical proxy ─────────────────────────────────────────────────────

def run_part_a(M, cl):
    """diagnosis_subtype — Non-Papillary (y=1) vs Papillary (y=0).

    Non-Papillary ≈ muscle-invasive / Basal-squamous molecular subtype
    Papillary     ≈ non-muscle-invasive / Luminal-papillary molecular subtype
    """
    print("\n── Part A: clinical proxy (Non-Papillary vs Papillary) ──", flush=True)
    labels = cl["diagnosis_subtype"].reindex(M.columns)
    keep   = labels.isin(["Non-Papillary", "Papillary"])
    Msub   = M.loc[:, keep].dropna(how="any")
    y      = (labels[keep] == "Non-Papillary").astype(int).values

    return _run(
        Msub, y,
        dataset_label  = "TCGA_BLCA_NonPapillary_vs_Papillary",
        endpoint_label = "clinical_subtype_nonpapillary_vs_papillary",
        expected_pole  = "?",   # fill after first run
        out_csv        = os.path.join(REPO, "external_validation_blca_clinical.csv"),
        novel_csv      = os.path.join(REPO, "novel_genes_blca_clinical.csv"),
    )


# ── Part B: Robertson 2017 molecular subtypes ──────────────────────────────────

def run_part_b(M, subtype_tsv):
    """Basal-squamous (y=1) vs Luminal-papillary (y=0) from Robertson et al. 2017 Cell."""
    print("\n── Part B: molecular subtypes (Robertson 2017) ──", flush=True)

    st = pd.read_csv(subtype_tsv, sep="\t")
    # normalise column names: accept 'sample'/'Sample'/'sampleID', 'subtype'/'Subtype'/'molecular_subtype'
    st.columns = [c.strip().lower().replace(" ", "_") for c in st.columns]
    sample_col  = next(c for c in st.columns if "sample" in c)
    subtype_col = next(c for c in st.columns if "subtype" in c)
    st = st.rename(columns={sample_col: "sample", subtype_col: "subtype"}).set_index("sample")

    print(f"  all subtypes in file: {dict(st['subtype'].value_counts())}")

    keep_map = {"Basal-squamous": 1, "Luminal-papillary": 0}
    mask     = st["subtype"].isin(keep_map)
    common   = st[mask].index.intersection(M.columns)
    Msub     = M.loc[:, common].dropna(how="any")
    y        = st.loc[common, "subtype"].map(keep_map).values

    return _run(
        Msub, y,
        dataset_label  = "TCGA_BLCA_BasalSquamous_vs_LuminalPapillary_Robertson2017",
        endpoint_label = "molecular_subtype_basal_vs_luminal",
        expected_pole  = "basal-squamous (KRT5/TP63/KRT14 expected — see lung overlap)",
        out_csv        = os.path.join(REPO, "external_validation_blca_subtype.csv"),
        novel_csv      = os.path.join(REPO, "novel_genes_blca_subtype.csv"),
    )


# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    cache_dir   = os.environ.get("BLCA_DIR", os.getcwd())
    subtype_tsv = os.environ.get("BLCA_SUBTYPE_TSV", "")

    print("== BLCA cross-cancer validation #4 ==")
    print(f"   cache dir:   {cache_dir}")
    if subtype_tsv:
        print(f"   subtype TSV: {subtype_tsv}")

    M, cl = _load_xena(cache_dir)
    print(f"   expression: {M.shape[0]} genes × {M.shape[1]} samples")

    results = {}
    results["a"] = run_part_a(M, cl)

    if subtype_tsv:
        if not os.path.exists(subtype_tsv):
            print(f"\n  [WARN] BLCA_SUBTYPE_TSV not found — skipping Part B")
        else:
            results["b"] = run_part_b(M, subtype_tsv)

    # ── cross-cancer summary ──────────────────────────────────────────────────
    print("\n" + "=" * 62)
    print("CROSS-CANCER SERIES SUMMARY")
    print("=" * 62)
    print("  #1 Lung  LUAD vs LUSC     : overlap 10/30, p=3.2e-16  → squamous pole")
    print("  #2 HNSC  tissue-indep     : panel AUROC 0.962")
    print("  #3 ESCA  ESCC vs EAC      : overlap  0/30, p=1.0      → adeno pole")
    if "a" in results:
        r = results["a"]
        print(f"  #4 BLCA  Non-Pap/Pap  (A): overlap {r['overlap']}/30, "
              f"p={r['p_hyper']:.2e}  | panel AUROC={r['basal_panel_auroc']}")
    if "b" in results:
        r = results["b"]
        print(f"  #4 BLCA  Basal/Luminal(B): overlap {r['overlap']}/30, "
              f"p={r['p_hyper']:.2e}  | panel AUROC={r['basal_panel_auroc']}")
    print()
    print("  INTERPRETATION GUIDE:")
    print("    overlap >> 0 + panel AUROC high → pan-epithelial basal axis confirmed")
    print("    overlap ≈  0 + panel AUROC high → panel transfers, axis present, but residual")
    print("                                       names a different pole (like ESCA adeno)")
    print("    overlap ≈  0 + panel AUROC low  → keratinization axis absent in this tissue")


if __name__ == "__main__":
    main()
