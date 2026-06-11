#!/usr/bin/env python3
"""Phase 1 — cross-study harmonization + meta-analysis on a real cohort.
Cohort: GSE57577 (engineered Dnmt3a2 in TKO mESC) + GSE77003 (Dnmt3a2 rescue in TKO mESC).
Two independent labs, different units (FPKM vs RPKM) and pipelines. We harmonize onto a
common gene-symbol space, correct study batch effects, and test whether the shared
biological axis — restoring Dnmt3a2 to TKO ESCs — reproduces across studies."""
import os, sys, tarfile, glob, gzip, shutil, numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from omniomics import geo, loaders, harmonize as hz, expression as ex
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.dirname(__file__); DATA = os.path.join(ROOT, "data"); FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

# ---- GSE57577: genotype-mean FPKM -> symbol matrix ----
d1 = os.path.join(DATA, "GSE57577")
tar = geo.download("GSE57575", "GSE57575_RAW.tar", d1)
with tarfile.open(tar) as t: t.extractall(d1)
for gz in glob.glob(os.path.join(d1, "*.gz")):
    if not os.path.exists(gz[:-3]):
        with gzip.open(gz) as fi, open(gz[:-3], "wb") as fo: shutil.copyfileobj(fi, fo)
mat, names, _ = loaders.load_cufflinks_fpkm_dir(d1)
genos = ["WT","WWD","R","TKO"]
gmean = pd.DataFrame({g: mat[[c for c in mat.columns if c.startswith(g+"_")]].mean(axis=1) for g in genos})
A = loaders.gene_symbol_index(gmean, names)          # symbol x {WT,WWD,R,TKO}
A.columns = [f"57577_{c}" for c in A.columns]

# ---- GSE77003: RPKM matrix -> symbol matrix ----
d2 = os.path.join(DATA, "GSE77003")
f2 = geo.download("GSE77003", "GSE77003_rpkm.txt.gz", d2)
B = loaders.load_matrix(f2, gene_col=0, drop_cols=["RefSeq"])   # index=Gene symbol
B = loaders.gene_symbol_index(B)
B.columns = [f"77003_{c.replace('mRNA_','')}" for c in B.columns]
print(f"[cohort] GSE57577 {A.shape} | GSE77003 {B.shape}")

# ---- harmonize onto common gene panel ----
studies = {"GSE57577": A, "GSE77003": B}
genes = hz.common_gene_panel(studies)
X, study_of = hz.assemble(studies, genes=genes, log=True)
Xqn = hz.quantile_normalize(X)
print(f"[cohort] common gene-symbol panel: {len(genes)} genes; combined {X.shape[1]} samples")

# batch effect before/after correction
bv_before, _ = hz.batch_variance_explained(Xqn, study_of)
Xc = hz.combat_lite(Xqn, study_of)
bv_after, _ = hz.batch_variance_explained(Xc, study_of)
print(f"[cohort] PC1 study-variance  before={bv_before['PC1']:.2f}  after={bv_after['PC1']:.2f}")
print(f"[cohort] PC2 study-variance  before={bv_before['PC2']:.2f}  after={bv_after['PC2']:.2f}")

# ---- shared axis: log2FC( Dnmt3a2-rescued / TKO ) per study ----
fc1 = np.log2(gmean["WT"]+1) - np.log2(gmean["TKO"]+1)        # GSE57577 WT(=Dnmt3a2) vs TKO
fc1 = pd.Series(fc1.values, index=[str(names.get(i,i)).upper() for i in gmean.index]).groupby(level=0).mean()
B_raw = loaders.gene_symbol_index(loaders.load_matrix(f2, gene_col=0, drop_cols=["RefSeq"]))
resc = B_raw[[c for c in B_raw.columns if "TKO3a2" in c]].mean(axis=1)
tko  = B_raw[[c for c in B_raw.columns if c=="mRNA_TKO"]].mean(axis=1)
fc2 = np.log2(resc+1) - np.log2(tko+1)                       # GSE77003 TKO3a2 vs TKO
mfc = pd.concat({"GSE57577": fc1, "GSE77003": fc2}, axis=1).dropna()
expr_ok = (B_raw[[c for c in B_raw.columns if "TKO" in c]].max(axis=1).reindex(mfc.index) > 1)
mfc = mfc[expr_ok.fillna(False)]
rho = stats.spearmanr(mfc["GSE57577"], mfc["GSE77003"]).correlation
r = stats.pearsonr(mfc["GSE57577"], mfc["GSE77003"])[0]
# directional concordance among genes strongly changed in study 1
strong = mfc[mfc["GSE57577"].abs() >= 1]
conc = (np.sign(strong["GSE57577"]) == np.sign(strong["GSE77003"])).mean()
print(f"[meta] shared Dnmt3a2-rescue axis on {len(mfc)} genes: Spearman={rho:.3f}, Pearson={r:.3f}")
print(f"[meta] sign concordance among |FC1|>=1 genes (n={len(strong)}): {100*conc:.1f}%")
mfc.to_csv(os.path.join(ROOT, "cross_study_rescue_FC.csv"))

# ---- figures ----
def scatter_pca(ax, M, title):
    pcs,_ = hz.pca(M, n=2)
    for sid,c in zip(study_of.unique(), ["#1f77b4","#d62728"]):
        m = (study_of.values==sid)
        ax.scatter(pcs[m,0], pcs[m,1], c=c, label=sid, s=45, edgecolor="k", linewidth=.4)
    ax.set_title(title); ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.legend(fontsize=7)
fig,axes = plt.subplots(1,2, figsize=(11,4.6))
scatter_pca(axes[0], Xqn, f"Before batch correction\nPC1 study-var={bv_before['PC1']:.2f}")
scatter_pca(axes[1], Xc,  f"After ComBat-lite\nPC1 study-var={bv_after['PC1']:.2f}")
plt.tight_layout(); plt.savefig(os.path.join(FIG,"cohort_batch_correction.png"),dpi=140); plt.close()

fig,ax = plt.subplots(figsize=(5.5,5.2))
ax.scatter(mfc["GSE57577"], mfc["GSE77003"], s=6, alpha=0.3, color="gray")
ax.axhline(0,lw=.5,c="k"); ax.axvline(0,lw=.5,c="k")
ax.set_xlabel("GSE57577  log2FC (WT-Dnmt3a2 / TKO)")
ax.set_ylabel("GSE77003  log2FC (TKO+Dnmt3a2 / TKO)")
ax.set_title(f"Cross-study reproducibility of Dnmt3a2 rescue\nSpearman={rho:.2f}, sign-concordance={100*conc:.0f}%")
plt.tight_layout(); plt.savefig(os.path.join(FIG,"cohort_cross_study_FC.png"),dpi=140); plt.close()
print("\n[cohort] saved figures + cross_study_rescue_FC.csv")
