#!/usr/bin/env python3
"""Real N>2 batch correction: 3 heterogeneous mouse studies (FPKM / RPKM / log-UQ counts),
harmonized and corrected with combat_lite vs EB-ComBat. Shows study batch removal on genuine
(not synthetic) multi-study data — the N>2 regime where correction is non-trivial."""
import os, sys, glob, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from omniomics import loaders, harmonize as hz, geo
import tarfile, gzip, shutil
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT=os.path.dirname(__file__); DATA=os.path.join(ROOT,"data"); FIG=os.path.join(ROOT,"figures")

# GSE57577 (FPKM, per sample)
d1=os.path.join(DATA,"GSE57577")
if not glob.glob(os.path.join(d1,"*.fpkm_tracking")):
    tar=geo.download("GSE57575","GSE57575_RAW.tar",d1)
    with tarfile.open(tar) as t: t.extractall(d1)
    for gz in glob.glob(os.path.join(d1,"*.gz")):
        if not os.path.exists(gz[:-3]):
            with gzip.open(gz) as fi,open(gz[:-3],"wb") as fo: shutil.copyfileobj(fi,fo)
matA,namesA,_=loaders.load_cufflinks_fpkm_dir(d1); A=loaders.gene_symbol_index(matA,namesA)
A=np.log2(A+1)
# GSE77003 (RPKM)
B=loaders.gene_symbol_index(loaders.load_matrix(os.path.join(DATA,"GSE77003","GSE77003_rpkm.txt.gz"),gene_col=0,drop_cols=["RefSeq"]))
B=np.log2(B+1)
# GSE316549 (already log-UQ)
C=loaders.gene_symbol_index(loaders.load_matrix(os.path.join(DATA,"GSE316549","GSE316549_loguq.counts.txt.gz"),gene_col=0))

studies={"GSE57577":A,"GSE77003":B,"GSE316549":C}
genes=hz.common_gene_panel(studies)
X,study=hz.assemble(studies,genes=genes,log=False)
X=X.dropna(); Xqn=hz.quantile_normalize(X)
print(f"[n3] common genes {len(genes)} ; samples per study: "+", ".join(f"{s}={ (study==s).sum() }" for s in study.unique()))

bv0,_=hz.batch_variance_explained(Xqn,study)
acc0=hz.batch_pred_accuracy(Xqn,study.values)
Xl=hz.combat_lite(Xqn,study); accL=hz.batch_pred_accuracy(Xl,study.values); bvL,_=hz.batch_variance_explained(Xl,study)
Xe=hz.combat_eb(Xqn,study.values); accE=hz.batch_pred_accuracy(Xe,study.values); bvE,_=hz.batch_variance_explained(Xe,study)
res=pd.DataFrame({"PC1_study_var":[bv0['PC1'],bvL['PC1'],bvE['PC1']],
                  "batch_recoverability":[acc0,accL,accE]},
                 index=["uncorrected","ComBat-lite","EB-ComBat"])
res.to_csv(os.path.join(ROOT,"mouse_n3_combat.csv"))
print("\n[verdict] 3-study mouse harmonization (chance batch acc = 0.33):")
print(res.round(3).to_string())

# PCA before/after
fig,ax=plt.subplots(1,2,figsize=(11,4.4)); pal={"GSE57577":"#1f77b4","GSE77003":"#d62728","GSE316549":"#2ca02c"}
for M,a,t in [(Xqn,ax[0],"Before correction"),(Xe,ax[1],"After EB-ComBat")]:
    pcs,_=hz.pca(M,n=2)
    for s,c in pal.items():
        m=(study.values==s); a.scatter(pcs[m,0],pcs[m,1],s=40,c=c,label=s,edgecolor="k",linewidth=.3)
    a.set_title(t); a.set_xlabel("PC1"); a.set_ylabel("PC2"); a.legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(FIG,"mouse_n3_combat.png"),dpi=140); plt.close()
print("[n3] saved figure + metrics")
