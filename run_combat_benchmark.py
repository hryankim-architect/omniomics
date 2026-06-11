#!/usr/bin/env python3
"""N>2 batch-correction benchmark: why covariate-aware EB-ComBat beats naive location/scale
when batch is CONFOUNDED with biology. Controlled ground truth on TCGA-BRCA LumA-vs-LumB
(a subtle axis). Inject 4 confounded batches + location/scale; correct 4 ways; measure
batch removal WITHIN subtype (isolating technical batch from the confounded biology) and
biology preservation (subtype AUROC)."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from omniomics import loaders, harmonize as hz
from omniomics import config
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BR=config.brca_tcga_dir()
ROOT=os.path.dirname(__file__); FIG=os.path.join(ROOT,"figures"); rng=np.random.default_rng(0)

coh=pd.read_csv(os.path.join(BR,"cohort_v2.tsv"),sep="\t").set_index("sample_id")
lab=coh[coh["group"].isin(["LumA","LumB"])]["group"]
rna=loaders.gene_symbol_index(loaders.load_matrix(os.path.join(BR,"HiSeqV2.gz"),gene_col=0))
S=[s for s in lab.index if s in rna.columns]; lab=lab.reindex(S)
var=rna[S].var(axis=1).sort_values(ascending=False); X0=rna.loc[var.index[:600], S]
y=(lab=="LumB").astype(int).values

K=4
pB={0:[.50,.30,.15,.05], 1:[.05,.15,.30,.50]}        # strong confounding batch<->subtype
batch=np.array([rng.choice(K, p=pB[int(yi)]) for yi in y])
Xb=X0.values.copy()
for b in range(K):
    ix=np.where(batch==b)[0]; loc=rng.normal(0,0.7,Xb.shape[0]); scale=rng.uniform(0.6,1.6,Xb.shape[0])
    gm=Xb[:,ix].mean(axis=1,keepdims=True); Xb[:,ix]=(Xb[:,ix]-gm)*scale[:,None]+gm+loc[:,None]
Xb=pd.DataFrame(Xb,index=X0.index,columns=S)
print(f"[combat] LumA-vs-LumB {len(S)} samples (LumB={y.sum()}) x {X0.shape[0]} genes; confounded 4 batches")
print(pd.crosstab(batch,y).to_string())

mod=pd.get_dummies(lab).values.astype(float)
def subtype_auc(M): return cross_val_score(LogisticRegression(max_iter=3000,C=0.3),
    StandardScaler().fit_transform(M.T.values), y, cv=StratifiedKFold(5,shuffle=True,random_state=1),
    scoring="roc_auc").mean()
def within_subtype_batch_acc(M):
    accs=[]
    for cls in (0,1):
        ix=np.where(y==cls)[0]
        accs.append(cross_val_score(RandomForestClassifier(60,random_state=0),
                    M.iloc[:,ix].T.values, batch[ix], cv=3).mean())
    return np.mean(accs)

methods={
 "uncorrected":            Xb,
 "ComBat-lite (no cov)":   hz.combat_lite(Xb, pd.Series(batch,index=S)),
 "EB-ComBat (no cov)":     hz.combat_eb(Xb, batch, mod=None),
 "EB-ComBat (+covariate)": hz.combat_eb(Xb, batch, mod=mod),
}
rows=[{"method":n,"within_subtype_batch_acc":within_subtype_batch_acc(M),"subtype_AUROC":subtype_auc(M)}
      for n,M in methods.items()]
res=pd.DataFrame(rows).set_index("method"); res.to_csv(os.path.join(ROOT,"combat_benchmark.csv"))
print("\n[verdict] N=4 confounded-batch correction:")
print(res.round(3).to_string())

fig,ax=plt.subplots(1,2,figsize=(12,4.4)); C=["#c44","#bb8844","#7fa8c9","#2E5A87"]
ax[0].bar(range(4),res["within_subtype_batch_acc"],color=C); ax[0].set_xticks(range(4)); ax[0].set_xticklabels(res.index,rotation=15,ha="right",fontsize=8)
ax[0].set_ylabel("within-subtype batch recoverability\n(lower = batch removed)"); ax[0].set_title("Technical batch removal")
ax[1].bar(range(4),res["subtype_AUROC"],color=C); ax[1].set_xticks(range(4)); ax[1].set_xticklabels(res.index,rotation=15,ha="right",fontsize=8)
ax[1].set_ylim(0.5,1.0); ax[1].set_ylabel("subtype AUROC (higher = biology kept)"); ax[1].set_title("Biology preserved under confounding")
for i,v in enumerate(res["subtype_AUROC"]): ax[1].text(i,v+0.01,f"{v:.3f}",ha="center",fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(FIG,"combat_benchmark.png"),dpi=140); plt.close()
print("[combat] saved figure + metrics")
