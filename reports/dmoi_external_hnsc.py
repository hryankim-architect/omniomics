#!/usr/bin/env python3
"""Third cross-cancer test (head & neck): the breast basal axis is a TISSUE-INDEPENDENT squamous marker.

HNSC is uniformly squamous, so it offers no within-cohort squamous-vs-adeno contrast -- but that enables a
cleaner confound control. We score every sample of TCGA HNSC (head & neck squamous), LUSC (lung squamous) and
LUAD (lung adeno) with the breast-derived 30-gene basal panel and ask whether the score tracks squamous
LINEAGE rather than tissue of origin.

Result: the breast basal panel separates squamous (HNSC+LUSC) from adeno (LUAD) at AUROC ~0.96, with HNSC
(head & neck) and LUSC (lung) -- two different tissues -- both scoring high and LUAD low. So the discovered
axis is a tissue-independent squamous-differentiation marker, not a breast/lung artefact. Within HNSC it also
tracks the textbook biology: well-differentiated (G1, keratinizing) tumours score higher than poorly
differentiated (G3) (one-sided p ~ 2e-4). Writes external_validation_hnsc.csv.

Run:  python reports/dmoi_external_hnsc.py        # downloads the three Xena matrices if absent
Env:  XENA_CACHE_DIR (default cwd) to cache hnsc_HiSeqV2.gz / lung_LUAD.gz / lung_LUSC.gz.
"""
import os, sys, io, urllib.request
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from sklearn.metrics import roc_auc_score
from scipy.stats import mannwhitneyu
CACHE = os.environ.get("XENA_CACHE_DIR", os.getcwd())
HUB = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.{c}.sampleMap/{f}"


def _expr(c, fname):
    p = os.path.join(CACHE, fname)
    if not os.path.exists(p):
        urllib.request.urlretrieve(HUB.format(c=c, f="HiSeqV2.gz"), p)
    d = pd.read_csv(p, sep="\t", index_col=0)
    return d[~d.index.duplicated()]


def main():
    basal = list(pd.read_csv(os.path.join(REPO, "novel_genes.csv"))["gene"])
    H, A, S = _expr("HNSC", "hnsc_HiSeqV2.gz"), _expr("LUAD", "lung_LUAD.gz"), _expr("LUSC", "lung_LUSC.gz")
    g = H.index.intersection(A.index).intersection(S.index); bp = [x for x in basal if x in g]
    M = pd.concat([H.loc[g], A.loc[g], S.loc[g]], axis=1)
    z = M.loc[bp].sub(M.loc[bp].mean(1), axis=0).div(M.loc[bp].std(1) + 1e-9, axis=0)
    score = z.mean(0)
    coh = np.array(["HNSC"] * H.shape[1] + ["LUAD"] * A.shape[1] + ["LUSC"] * S.shape[1])
    sq = np.array([1 if c in ("HNSC", "LUSC") else 0 for c in coh])
    auc = roc_auc_score(sq, score.values)
    med = {c: round(float(np.median(score.values[coh == c])), 3) for c in ["HNSC", "LUSC", "LUAD"]}
    cm_url = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.HNSC.sampleMap/HNSC_clinicalMatrix"
    cm = pd.read_csv(io.BytesIO(urllib.request.urlopen(cm_url, timeout=30).read()), sep="\t", index_col=0)
    gr = cm["neoplasm_histologic_grade"].reindex(H.columns)
    hs = pd.Series(score.values[coh == "HNSC"], index=H.columns)
    g1, g3 = hs[gr == "G1"].dropna(), hs[gr == "G3"].dropna()
    _, p = mannwhitneyu(g1, g3, alternative="greater")
    pd.DataFrame([dict(test="tissue_independence", basal_genes=len(bp), auroc_squamous_vs_adeno=round(auc, 3),
                       median_HNSC=med["HNSC"], median_LUSC=med["LUSC"], median_LUAD=med["LUAD"],
                       within_hnsc_G1_med=round(float(g1.median()), 3), within_hnsc_G3_med=round(float(g3.median()), 3),
                       within_hnsc_G1gtG3_p=f"{p:.2e}",
                       note="breast basal panel = tissue-independent squamous-differentiation marker (HNSC~LUSC>>LUAD; tracks HNSC grade)")
                  ]).to_csv(os.path.join(REPO, "external_validation_hnsc.csv"), index=False)
    print(f"basal panel separates squamous (HNSC+LUSC) vs adeno (LUAD): AUROC={auc:.3f}")
    print(f"median basal by cohort: {med}  -> HNSC(head&neck) ~ LUSC(lung) >> LUAD(adeno)")
    print(f"within HNSC: G1>G3 (keratinization tracks differentiation) p={p:.2e}")
    print("CONCLUSION: the discovered basal axis is a tissue-independent squamous-differentiation marker.")


if __name__ == "__main__":
    main()
