#!/usr/bin/env python3
"""Independent validation of the dmoi-brca-poc thesis on the LumA-vs-LumB axis.
Hypothesis: where the discriminating biology is proliferation/ER (LumA vs LumB) and priors are
POLE-SPECIFIC, a compact multi-omics representation should rival/beat a big RNA model — unlike
HER2-vs-Luminal where naive fusion failed. Pole priors (from the project's own pathway.py):
  LumA pole = ESTROGEN_RESPONSE_EARLY/LATE ; LumB pole = E2F_TARGETS/G2M_CHECKPOINT/MYC_TARGETS_V1."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from omniomics import loaders, methylation as me, multiomics as mo
from omniomics import config
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold

BR=config.brca_tcga_dir()
MT=config.cache_dir(); GMT=config.hallmark_gmt()
ROOT=os.path.dirname(__file__); FIG=os.path.join(ROOT,"figures")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# labels: LumA vs LumB (cohort_v2), RNA available
coh=pd.read_csv(os.path.join(BR,"cohort_v2.tsv"),sep="\t").set_index("sample_id")
lab=coh[coh["group"].isin(["LumA","LumB"])]["group"]

rna=loaders.gene_symbol_index(loaders.load_matrix(os.path.join(BR,"HiSeqV2.gz"),gene_col=0))
meth=me.load_hm450_promoter(os.path.join(MT,"promoter_meth.tsv.gz"),os.path.join(MT,"promoter_probe_gene.tsv"))
S=[s for s in lab.index if s in rna.columns and s in meth.columns]
G=sorted(set(rna.index)&set(meth.index))
rna=rna.loc[G,S]; meth=meth.loc[G,S]; y=(lab.reindex(S)=="LumB").astype(int).values
print(f"[lumab] {len(S)} samples (LumA={ (y==0).sum() }, LumB={ (y==1).sum() }) x {len(G)} genes")

gmt=mo.load_gmt(GMT)
POLE={"LumA":["HALLMARK_ESTROGEN_RESPONSE_EARLY","HALLMARK_ESTROGEN_RESPONSE_LATE"],
      "LumB":["HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MYC_TARGETS_V1"]}
def pole_genes(p):
    g=set()
    for s in POLE[p]: g|=set(gmt.get(s,[]))
    return [x for x in g if x in rna.index]
def zscore(df): return pd.DataFrame(StandardScaler().fit_transform(df.T.fillna(df.T.mean()).values),
                                    index=df.columns, columns=df.index)
zR=zscore(rna); zM=zscore(-meth)   # -meth: high methylation -> low expression direction

def pole_feats(zmat):
    return np.column_stack([zmat[pole_genes(p)].mean(axis=1).values for p in ["LumA","LumB"]])
RNA_pole=pole_feats(zR); METH_pole=pole_feats(zM)
DIS=RNA_pole-METH_pole                                   # per-pole RNA-vs-meth disagreement (the DMOI scalar idea)

cv=StratifiedKFold(5,shuffle=True,random_state=0)
def auc(X): return cross_val_score(LogisticRegression(max_iter=4000,C=1.0),X,y,cv=cv,scoring="roc_auc").mean()
var=rna.var(axis=1).sort_values(ascending=False); top=rna.loc[var.index[:2000]]
Xr=StandardScaler().fit_transform(top.T.values)

res={
 "RNA (2000 genes)":                 auc(Xr),
 "RNA pole (2 feats)":               auc(StandardScaler().fit_transform(RNA_pole)),
 "METH pole (2 feats)":              auc(StandardScaler().fit_transform(METH_pole)),
 "RNA+METH pole (4)":                auc(StandardScaler().fit_transform(np.hstack([RNA_pole,METH_pole]))),
 "RNA+METH pole + disagree (6)":     auc(StandardScaler().fit_transform(np.hstack([RNA_pole,METH_pole,DIS]))),
}
print("\n[verdict] LumA-vs-LumB AUROC (5-fold CV):")
for k,v in res.items(): print(f"   {k:34s} {v:.3f}")
base=res["RNA pole (2 feats)"]; best_mo=max(res["RNA+METH pole (4)"],res["RNA+METH pole + disagree (6)"])
print(f"\n[gain] multi-omics pole vs RNA pole: {best_mo-base:+.3f}  "
      f"({'GAIN' if best_mo>base else 'no gain'})")
pd.Series(res).to_csv(os.path.join(ROOT,"lumab_dmoi_auroc.csv"))

fig,ax=plt.subplots(figsize=(8,3.8))
ks=list(res); vs=[res[k] for k in ks]
ax.barh(range(len(ks)),vs,color=["#999","#88a","#bb8844","#4488aa","#2E5A87"])
ax.set_yticks(range(len(ks))); ax.set_yticklabels(ks,fontsize=8); ax.set_xlim(0.5,1.0)
ax.axvline(res["RNA pole (2 feats)"],ls="--",c="k",lw=.7,label="RNA-pole baseline")
for i,v in enumerate(vs): ax.text(v+0.003,i,f"{v:.3f}",va="center",fontsize=8)
ax.set_xlabel("LumA-vs-LumB AUROC (5-fold CV)")
ax.set_title("Pole-conditioned multi-omics on the LumA/LumB axis"); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(FIG,"lumab_dmoi_auroc.png"),dpi=140); plt.close()
print("[lumab] saved figure + metrics")
