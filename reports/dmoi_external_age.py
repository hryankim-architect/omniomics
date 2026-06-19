#!/usr/bin/env python3
"""External validation of the anchored integrator on the textbook methylation endpoint
(epigenetic-clock AGE) in TCGA-BRCA normal-adjacent tissue (sample type -11, no tumour-
proliferation confound). Honest result: even here a 1,500-gene RNA model beats a genome-wide
3,000-CpG methylation model (and a CpG-selected, clock-mimicking one), so the anchored gate
correctly stays at the RNA anchor -- genuine multi-omics gains need *curated* methylation
biomarkers (the real clock uses 353 specific CpGs). The frame still adapts: anchoring on
methylation, RNA earns a small gain. Demonstrates `anchored_integrate` on real data.

Inputs (built earlier in the session, or regenerate): rna1500.tsv (1500-gene panel),
meth_gw.tsv (~3000 genome-wide CpGs x all 450K samples), TCGA dir with HiSeqV2.gz,
BRCA_clinicalMatrix.tsv. Set BRCA_DIR / OUT env or edit below.
"""
import os, sys, gzip, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from omniomics import multiomics as mo
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score

D = os.environ.get("BRCA_DIR", "/path/to/tcga_brca")
OUT = os.environ.get("OUT", os.path.dirname(os.path.abspath(__file__)))

mgw = pd.read_csv(f"{OUT}/meth_gw.tsv", sep="\t", index_col=0)
panel = list(pd.read_csv(f"{OUT}/rna1500.tsv", sep="\t", index_col=0).index)
with gzip.open(f"{D}/HiSeqV2.gz", "rt") as fh:
    rcols = fh.readline().rstrip("\n").split("\t")[1:]
cl = pd.read_csv(f"{D}/BRCA_clinicalMatrix.tsv", sep="\t",
                 usecols=["sampleID", "age_at_initial_pathologic_diagnosis"])
cl["age"] = pd.to_numeric(cl["age_at_initial_pathologic_diagnosis"], errors="coerce")
age = dict(zip(cl.sampleID, cl.age))

norm = [s for s in mgw.columns if s.endswith("-11") and s in set(rcols) and age.get(s) == age.get(s)]
rna = pd.read_csv(f"{D}/HiSeqV2.gz", sep="\t", index_col=0, usecols=["sample"] + norm)
rna = rna.loc[[g for g in panel if g in rna.index]]
S = [s for s in norm if s in rna.columns]
y = (np.array([age[s] for s in S]) > np.median([age[s] for s in S])).astype(int)
Xr = rna[S].T.values
Mg = mgw[S].T.astype(float); Mg = Mg.fillna(Mg.mean()).fillna(0.0); Xm = Mg.values

cv = StratifiedKFold(5, shuffle=True, random_state=0)
lr = lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000, C=0.3))
sel = lambda: make_pipeline(SelectKBest(f_classif, k=120), StandardScaler(),
                            LogisticRegression(max_iter=5000, C=0.3))
oof = lambda est, X: cross_val_predict(est, X, y, cv=cv, method="predict_proba")[:, 1]

print(f"TCGA-BRCA normal-adjacent AGE (older-than-median)  n={len(S)}  RNA genes={Xr.shape[1]}  CpGs={Xm.shape[1]}")
print(f"  RNA alone                      AUROC={roc_auc_score(y, oof(lr(), Xr)):.3f}")
print(f"  Methylation (genome-wide)      AUROC={roc_auc_score(y, oof(lr(), Xm)):.3f}")
print(f"  Methylation (CpG-selected/clock) AUROC={roc_auc_score(y, oof(sel(), Xm)):.3f}")
dd = [mo.anchored_integrate(Xr, Xm, y, cv=5, random_state=s)["delta"] for s in range(6)]
print(f"  RNA-anchored + meth (anchored_integrate, 6 seeds): mean delta={np.mean(dd):+.3f}")
d2 = [mo.anchored_integrate(Xm, Xr, y, cv=5, random_state=s)["delta"] for s in range(6)]
print(f"  METH-anchored + RNA (frame adapts):                mean delta={np.mean(d2):+.3f}")
