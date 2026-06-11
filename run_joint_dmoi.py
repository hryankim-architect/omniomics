#!/usr/bin/env python3
"""(1) MOFA-style joint RNA+methylation embedding with per-view variance decomposition.
(2) DMOI-lite: pathway-conditioned fusion to test whether STRUCTURED multi-omics beats
    RNA-alone, where naive concatenation did not."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from omniomics import loaders, methylation as me, multiomics as mo
from omniomics import config
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BR=config.brca_tcga_dir()
MT=config.cache_dir(); GMT=config.hallmark_gmt()
ROOT=os.path.dirname(__file__); FIG=os.path.join(ROOT,"figures")

rna = loaders.gene_symbol_index(loaders.load_matrix(os.path.join(BR,"HiSeqV2.gz"), gene_col=0))
meth = me.load_hm450_promoter(os.path.join(MT,"promoter_meth.tsv.gz"), os.path.join(MT,"promoter_probe_gene.tsv"))
S = sorted(set(rna.columns)&set(meth.columns)); G = sorted(set(rna.index)&set(meth.index))
rna=rna.loc[G,S]; meth=meth.loc[G,S]
coh=pd.read_csv(os.path.join(BR,"cohort_v4.tsv"),sep="\t").set_index("sample_id")["group"].reindex(S)
keep=coh.isin(["HER2","Luminal"]); y=(coh[keep]=="HER2").astype(int).values
Sk=[s for s,k in zip(S,keep) if k]
cv=StratifiedKFold(5,shuffle=True,random_state=0); clf=lambda: LogisticRegression(max_iter=4000,C=0.1)
def auc(X): return cross_val_score(clf(),X,y,cv=cv,scoring="roc_auc").mean()
print(f"[joint] matched {len(S)} samples x {len(G)} genes")

# ===== (1) MOFA-style joint embedding =====
T, R2 = mo.mofa_lite({"RNA": rna, "Methylation": meth}, k=10)
R2.to_csv(os.path.join(ROOT,"mofa_variance_explained.csv"))
shared = R2.loc["RNA"].gt(0.002) & R2.loc["Methylation"].gt(0.002)
print("[MOFA] per-view variance explained (top factors):")
print((R2.T*100).round(2).to_string())
print(f"[MOFA] shared factors (both views >0.2%): {[f for f in R2.columns if shared[f]]}")
auc_joint = auc(StandardScaler().fit_transform(T.loc[Sk].values))
print(f"[MOFA] subtype AUROC from 10 joint factors only: {auc_joint:.3f}")

# ===== (2) DMOI-lite pathway-conditioned fusion =====
gmt = mo.load_gmt(GMT)
rna_pw  = mo.pathway_scores(rna,  gmt)          # pathways x samples
meth_pw = mo.pathway_scores(-meth, gmt)         # -meth aligns with expression direction
common_pw = sorted(set(rna_pw.index)&set(meth_pw.index))
RP = rna_pw.loc[common_pw, Sk].T.values         # samples x pathways
MP = meth_pw.loc[common_pw, Sk].T.values
DIS = StandardScaler().fit_transform(RP) - StandardScaler().fit_transform(MP)  # per-pathway RNA-vs-meth disagreement
print(f"[DMOI] {len(common_pw)} Hallmark pathways scored per omics")

# baselines (raw features)
var=rna[Sk].var(axis=1).sort_values(ascending=False); top=var.index[:2000]
Xr=StandardScaler().fit_transform(rna.loc[top,Sk].T.values)
Xm=StandardScaler().fit_transform(meth.loc[top,Sk].fillna(meth.loc[top,Sk].mean()).T.values)
res={
 "RNA (2000 genes)":            auc(Xr),
 "RNA+METH concat (naive)":     auc(np.hstack([Xr,Xm])),
 "RNA pathway (50)":            auc(StandardScaler().fit_transform(RP)),
 "RNA+METH pathway":            auc(StandardScaler().fit_transform(np.hstack([RP,MP]))),
 "RNA+METH pathway + disagree": auc(StandardScaler().fit_transform(np.hstack([RP,MP,DIS]))),
 "MOFA joint factors (10)":     auc_joint,
}
print("\n[verdict] subtype AUROC (5-fold CV):")
for k,v in res.items(): print(f"   {k:32s} {v:.3f}")
pd.Series(res).to_csv(os.path.join(ROOT,"dmoi_fusion_auroc.csv"))

# ===== figures =====
fig,ax=plt.subplots(figsize=(8,2.8))
im=ax.imshow((R2.values*100), aspect="auto", cmap="viridis")
ax.set_yticks([0,1]); ax.set_yticklabels(R2.index); ax.set_xticks(range(R2.shape[1])); ax.set_xticklabels(R2.columns)
for i in range(R2.shape[0]):
    for j in range(R2.shape[1]): ax.text(j,i,f"{R2.values[i,j]*100:.1f}",ha="center",va="center",color="w",fontsize=7)
plt.colorbar(im,label="% variance explained"); ax.set_title("MOFA-style joint factors: variance explained per omics")
plt.tight_layout(); plt.savefig(os.path.join(FIG,"mofa_variance.png"),dpi=140); plt.close()

fig,ax=plt.subplots(figsize=(8,4.2))
ks=list(res); vs=[res[k] for k in ks]
colors=["#999","#c44","#88a","#4488aa","#2E5A87","#6a3d9a"]
ax.barh(range(len(ks)), vs, color=colors); ax.set_yticks(range(len(ks))); ax.set_yticklabels(ks,fontsize=8)
ax.set_xlim(0.5,1.0); ax.axvline(res["RNA (2000 genes)"],ls="--",c="k",lw=.7,label="RNA-only baseline")
for i,v in enumerate(vs): ax.text(v+0.003,i,f"{v:.3f}",va="center",fontsize=8)
ax.set_xlabel("subtype AUROC (5-fold CV)"); ax.set_title("Does structured multi-omics beat RNA alone?"); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(FIG,"dmoi_fusion_auroc.png"),dpi=140); plt.close()
print("\n[joint] saved figures + metrics")
