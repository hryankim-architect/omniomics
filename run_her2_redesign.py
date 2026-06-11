#!/usr/bin/env python3
"""HER2 amplicon prior redesign. Phase 6b found methylation HURTS HER2-vs-Luminal because the
HER2 pole is the ERBB2 amplicon (copy-number/transcriptional, not methylation-regulated).
Redesign: route methylation SELECTIVELY — only to the ER/luminal pole (where it's informative),
RNA-only on the amplicon pole. Test whether selective routing avoids the hurt."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from omniomics import loaders, methylation as me, multiomics as mo
from omniomics import config
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BR=config.brca_tcga_dir(); MT=config.cache_dir()
ROOT=os.path.dirname(__file__); FIG=os.path.join(ROOT,"figures")

rna=loaders.gene_symbol_index(loaders.load_matrix(os.path.join(BR,"HiSeqV2.gz"),gene_col=0))
Mm=pd.read_csv(os.path.join(MT,"enh_meth_all.tsv.gz"),sep="\t",index_col=0).astype("float32")
meth=me.aggregate(Mm, os.path.join(MT,"pg_enhancer.tsv")); del Mm
coh=pd.read_csv(os.path.join(BR,"cohort_v4.tsv"),sep="\t").set_index("sample_id")
lab=coh[coh["group"].isin(["HER2","Luminal"])]["group"]
S=[s for s in lab.index if s in rna.columns and s in meth.columns]
y=(lab.reindex(S)=="HER2").astype(int).values
G=sorted(set(rna.index)&set(meth.index)); R=rna.loc[G,S]; Me=meth.loc[G,S]
MARK={"HER2":["ERBB2","GRB7","STARD3","PGAP3","MIEN1"],"ER":["ESR1","GATA3","FOXA1","XBP1","PGR"]}
D=mo.dmoi_representation(R,Me,{"HER2":["HER2"],"ER":["ER"]},MARK)
print(f"[her2] {len(S)} samples (HER2={y.sum()})")

feats={
 "RNA pole (2)":                D[["rna_HER2","rna_ER"]].values,
 "naive DMOI (meth both, 4)":   D[["rna_HER2","rna_ER","meth_HER2","meth_ER"]].values,
 "SELECTIVE (meth ER only, 3)": D[["rna_HER2","rna_ER","meth_ER"]].values,
 "SELECTIVE + disagree_ER (4)": D[["rna_HER2","rna_ER","meth_ER","disagree_ER"]].values,
}
def rcv(X,y,rep=15):
    a=[]
    for r in range(rep):
        for tr,te in StratifiedKFold(5,shuffle=True,random_state=r).split(np.zeros(len(y)),y):
            sc=StandardScaler().fit(X[tr]); clf=LogisticRegression(max_iter=4000,C=1.0).fit(sc.transform(X[tr]),y[tr])
            a.append(roc_auc_score(y[te],clf.predict_proba(sc.transform(X[te]))[:,1]))
    return np.array(a)
Rr={k:rcv(v,y) for k,v in feats.items()}
summary=pd.DataFrame({"AUROC_mean":{k:v.mean() for k,v in Rr.items()},
                      "AUROC_std":{k:v.std() for k,v in Rr.items()}})
summary.to_csv(os.path.join(ROOT,"her2_redesign_auroc.csv"))
print("\n[verdict] HER2-vs-Luminal, 15x5 CV:")
print(summary.round(4).to_string())
base="RNA pole (2)"
for k in feats:
    if k==base: continue
    d=Rr[k].mean()-Rr[base].mean(); p=stats.wilcoxon(Rr[k],Rr[base]).pvalue
    print(f"  {k:30s} vs RNA-pole: {d:+.4f}  p={p:.2e}  ({'better' if d>0 and p<0.05 else 'worse' if d<0 and p<0.05 else 'ns'})")
d=Rr["SELECTIVE (meth ER only, 3)"].mean()-Rr["naive DMOI (meth both, 4)"].mean()
p=stats.wilcoxon(Rr["SELECTIVE (meth ER only, 3)"],Rr["naive DMOI (meth both, 4)"]).pvalue
print(f"  selective vs naive routing: {d:+.4f}  p={p:.2e}")

fig,ax=plt.subplots(figsize=(8.5,4)); ks=list(feats); m=[Rr[k].mean() for k in ks]; e=[Rr[k].std() for k in ks]
ax.bar(range(len(ks)),m,yerr=e,capsize=4,color=["#999","#c44","#2E5A87","#4488aa"])
ax.axhline(Rr[base].mean(),ls="--",c="k",lw=.7,label="RNA-pole baseline")
ax.set_xticks(range(len(ks))); ax.set_xticklabels(ks,rotation=12,ha="right",fontsize=8); ax.set_ylim(0.88,0.96)
for i,(mm,ee) in enumerate(zip(m,e)): ax.text(i,mm+ee+0.001,f"{mm:.3f}",ha="center",fontsize=8)
ax.set_ylabel("HER2-vs-Luminal AUROC (15x5 CV)"); ax.set_title("HER2 prior redesign: selective methylation routing"); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(FIG,"her2_redesign.png"),dpi=140); plt.close()
print("[her2] saved figure + metrics")
