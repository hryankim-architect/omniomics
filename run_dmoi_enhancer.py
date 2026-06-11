#!/usr/bin/env python3
"""Formally incorporate ENHANCER methylation into DMOI-style structured fusion and test —
rigorously — whether it yields a real multi-omics gain over RNA alone on LumA-vs-LumB.
DMOI representation = per-pole {RNA score, methylation score, RNA-vs-meth disagreement}.
Significance via repeated 5-fold CV (20 repeats, shared splits) + paired Wilcoxon."""
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

BR=config.brca_tcga_dir()
MT=config.cache_dir(); GMT=config.hallmark_gmt()
ROOT=os.path.dirname(__file__); FIG=os.path.join(ROOT,"figures")

coh=pd.read_csv(os.path.join(BR,"cohort_v2.tsv"),sep="\t").set_index("sample_id")
lab=coh[coh["group"].isin(["LumA","LumB"])]["group"]
rna=loaders.gene_symbol_index(loaders.load_matrix(os.path.join(BR,"HiSeqV2.gz"),gene_col=0))
M=pd.read_csv(os.path.join(MT,"ctx_meth_lumab.tsv.gz"),sep="\t",index_col=0).astype("float32")
meth_enh=me.aggregate(M, os.path.join(MT,"pg_enhancer.tsv"))
meth_pro=me.aggregate(M, os.path.join(MT,"pg_promoter.tsv"))
del M

S=[s for s in lab.index if s in rna.columns and s in meth_enh.columns and s in meth_pro.columns]
y=(lab.reindex(S)=="LumB").astype(int).values
rna=rna[S]; gmt=mo.load_gmt(GMT)
POLE={"LumA":["HALLMARK_ESTROGEN_RESPONSE_EARLY","HALLMARK_ESTROGEN_RESPONSE_LATE"],
      "LumB":["HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MYC_TARGETS_V1"]}
print(f"[dmoi-enh] {len(S)} samples (LumB={y.sum()}) ; building DMOI representations")

def align(meth):
    G=sorted(set(rna.index)&set(meth.index)); return rna.loc[G], meth.loc[G,S]
rE,mE=align(meth_enh); rP,mP=align(meth_pro)
D_enh=mo.dmoi_representation(rE,mE,POLE,gmt)      # samples x 6
D_pro=mo.dmoi_representation(rP,mP,POLE,gmt)

# feature matrices (all aligned to S, y)
feats={
 "RNA pole (2)":                 D_enh[[f"rna_{p}" for p in POLE]].values,
 "RNA+enhMeth pole (4)":         D_enh[[f"{a}_{p}" for a in ["rna","meth"] for p in POLE]].values,
 "DMOI enhancer (6: +disagree)": D_enh.values,
 "DMOI promoter (6, contrast)":  D_pro.values,
}

def repeated_cv(Xdict, y, repeats=20):
    out={k:[] for k in Xdict}
    for r in range(repeats):
        splits=list(StratifiedKFold(5,shuffle=True,random_state=r).split(np.zeros(len(y)),y))
        for k,X in Xdict.items():
            a=[]
            for tr,te in splits:
                sc=StandardScaler().fit(X[tr])
                clf=LogisticRegression(max_iter=5000,C=1.0).fit(sc.transform(X[tr]),y[tr])
                a.append(roc_auc_score(y[te], clf.predict_proba(sc.transform(X[te]))[:,1]))
            out[k].append(np.mean(a))
    return {k:np.array(v) for k,v in out.items()}

R=repeated_cv(feats,y,repeats=20)
# RNA 2000-gene reference (single pass, for context)
var=rna.var(axis=1).sort_values(ascending=False)
_sp=list(StratifiedKFold(5,shuffle=True,random_state=0).split(np.zeros(len(y)),y))
_X=rna.loc[var.index[:1500]].T.values
rna2000=np.mean([roc_auc_score(y[te], LogisticRegression(max_iter=2000,C=0.5).fit(
    StandardScaler().fit_transform(_X[tr]),y[tr]).predict_proba(
    StandardScaler().fit(_X[tr]).transform(_X[te]))[:,1]) for tr,te in _sp])
R["RNA ~1500 genes"]=np.array([rna2000])
summary=pd.DataFrame({"AUROC_mean":{k:v.mean() for k,v in R.items()},
                      "AUROC_std":{k:v.std() for k,v in R.items()}})
print("\n[verdict] LumA-vs-LumB, 20x5-fold repeated CV:")
print(summary.round(4).to_string())

base="RNA pole (2)"
for k in ["RNA+enhMeth pole (4)","DMOI enhancer (6: +disagree)","DMOI promoter (6, contrast)"]:
    d=R[k]-R[base]; p=stats.wilcoxon(R[k],R[base]).pvalue
    print(f"  {k:30s} gain={d.mean():+.4f}  paired Wilcoxon p={p:.2e}  ({'SIG' if p<0.05 and d.mean()>0 else 'ns/neg'})")
# enhancer vs promoter DMOI
d=R["DMOI enhancer (6: +disagree)"]-R["DMOI promoter (6, contrast)"]
print(f"  DMOI enhancer vs promoter: gain={d.mean():+.4f}  p={stats.wilcoxon(R['DMOI enhancer (6: +disagree)'],R['DMOI promoter (6, contrast)']).pvalue:.2e}")
summary.to_csv(os.path.join(ROOT,"dmoi_enhancer_auroc.csv"))

# figure
order=["RNA pole (2)","RNA+enhMeth pole (4)","DMOI enhancer (6: +disagree)","DMOI promoter (6, contrast)","RNA ~1500 genes"]
m=[R[k].mean() for k in order]; e=[R[k].std() for k in order]
fig,ax=plt.subplots(figsize=(8.5,4.2))
ax.bar(range(len(order)),m,yerr=e,capsize=4,
       color=["#999","#bb8844","#2E5A87","#7fa8c9","#555"])
ax.axhline(R["RNA pole (2)"].mean(),ls="--",c="k",lw=.7,label="RNA-pole baseline")
ax.set_xticks(range(len(order))); ax.set_xticklabels(order,rotation=15,ha="right",fontsize=8)
ax.set_ylim(0.85,0.96); ax.set_ylabel("LumA-vs-LumB AUROC (20x5 CV)")
ax.set_title("Enhancer methylation in DMOI-style fusion — does structure earn a gain?")
for i,(mm,ee) in enumerate(zip(m,e)): ax.text(i,mm+ee+0.001,f"{mm:.3f}",ha="center",fontsize=8)
ax.legend(fontsize=8); plt.tight_layout(); plt.savefig(os.path.join(FIG,"dmoi_enhancer.png"),dpi=140); plt.close()
print("[dmoi-enh] saved figure + metrics")
