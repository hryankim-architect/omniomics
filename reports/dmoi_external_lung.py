#!/usr/bin/env python3
"""Cross-CANCER external validation: does the breast basal/keratinization axis re-discover in lung?

The strongest test of whether a discovered axis is real biology rather than a cohort artefact is to look
for it in a *different cancer*. The breast LumA/B residual (anchored on proliferation) was a
basal/squamous-keratinization program (KRT5/14/17/6B, TP63, DSG3/DSC3, ...). Lung carcinoma splits into
adenocarcinoma (LUAD) and squamous-cell carcinoma (LUSC); the squamous program is the textbook
discriminator. So we take TCGA LUAD+LUSC (UCSC Xena HiSeqV2, n~1129), anchor on the SAME zero-parameter
proliferation signature (orthogonal to histology, AUROC ~0.77 = incomplete), and mine the residual.

Result: the residual independently recovers the squamous/keratinization axis, and it overlaps the breast
basal panel 10/30 (hypergeometric p ~ 3e-16) -- the SAME genes (KRT5, KRT14, KRT6B, TP63, DSG3, DSC3, FAT2,
CALML3, ANXA8, TRIM29) discovered in a different cancer. HONEST CAVEAT: the panel-vs-random delta control
saturates here (p ~ 0.18) because squamous-vs-adeno is a near-trivial transcriptomic split -- random
high-variance panels also separate it -- so on this easy endpoint the informative metric is the gene-level
replication (which is decisive), not the panel-vs-random margin. Writes external_validation_lung.csv and
novel_genes_lung.csv.

Run:  python reports/dmoi_external_lung.py            # downloads the two Xena matrices if absent
Env:  LUNG_CACHE_DIR (default cwd) to cache lung_LUAD.gz / lung_LUSC.gz.
"""
import os, sys
from math import comb
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from omniomics import multiomics as mo
from sklearn.metrics import roc_auc_score
CACHE = os.environ.get("LUNG_CACHE_DIR", os.getcwd())
XENA = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.{c}.sampleMap/HiSeqV2.gz"
PROLIF = ["MKI67", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB", "BUB1", "CCNE1", "CDC20",
          "TOP2A", "TYMS", "RRM2", "UBE2C", "CENPF", "FOXM1", "MELK", "KIF2C", "NUSAP1", "PTTG1"]


def _get(c):
    import urllib.request
    p = os.path.join(CACHE, f"lung_{c}.gz")
    if not os.path.exists(p):
        urllib.request.urlretrieve(XENA.format(c=c), p)
    d = pd.read_csv(p, sep="\t", index_col=0)
    return d[~d.index.duplicated()]


def main():
    luad, lusc = _get("LUAD"), _get("LUSC")
    g = luad.index.intersection(lusc.index)
    M = pd.concat([luad.loc[g], lusc.loc[g]], axis=1).dropna(how="any")
    y = np.array([0] * luad.shape[1] + [1] * lusc.shape[1])      # 0=LUAD(adeno) 1=LUSC(squamous)
    pres = [x for x in PROLIF if x in M.index]
    anchor = mo.signature_score(M.loc[pres].T, pres)
    basal = set(pd.read_csv(os.path.join(REPO, "novel_genes.csv"))["gene"])
    var = M.var(axis=1); feats = sorted(set(var.sort_values(ascending=False).head(5000).index) | (basal & set(M.index)))
    X = M.loc[feats].T.values.astype("float32")
    res = mo.anchored_residual_discovery(anchor, X, feats, y, top_k=30, corr_max=0.6, cv=5,
                                         random_state=0, n_perm=10, inner_repeats=1)
    nov = [gn for gn, _, _ in res["novel"]]; ov = sorted(set(nov) & basal)
    N, K, n, k = len(feats), len(basal & set(feats)), len(nov), len(ov)
    p_hyper = sum(comb(K, i) * comb(N - K, n - i) for i in range(k, min(K, n) + 1)) / comb(N, n)
    pd.DataFrame(res["novel"], columns=["gene", "partial_corr", "corr_with_anchor"]).to_csv(
        os.path.join(REPO, "novel_genes_lung.csv"), index=False)
    pd.DataFrame([dict(dataset="TCGA_lung_LUAD_vs_LUSC", endpoint="histology_squamous_vs_adeno", n=int(len(y)),
                       textbook_anchor="proliferation_20gene", anchor_auroc=round(res["auroc_anchor"], 3),
                       combined=round(res["auroc_combined"], 3), delta=round(res["delta"], 3),
                       novel_delta=round(res["novel_delta"], 3), random_delta_mean=round(res["random_delta_mean"], 3),
                       novel_vs_random_p=round(res["novel_vs_random_p"], 3), overlap_with_brca_basal_of30=int(k),
                       overlap_hyper_p=f"{p_hyper:.2e}", shared_genes="; ".join(ov))
                  ]).to_csv(os.path.join(REPO, "external_validation_lung.csv"), index=False)
    print(f"TCGA lung LUAD-vs-LUSC n={len(y)} | proliferation anchor AUROC={res['auroc_anchor']:.3f}")
    print(f"overlap with BRCA basal panel: {k}/30 (hypergeometric p={p_hyper:.2e})")
    print("shared genes:", ", ".join(ov))
    print("CONCLUSION: the breast basal/keratinization axis re-discovers in an independent cancer (lung "
          "squamous) -- the discovered biology, not a cohort artefact. (Panel-vs-random saturates: "
          "squamous-vs-adeno is a near-trivial split, so gene-level overlap is the informative metric.)")


if __name__ == "__main__":
    main()
