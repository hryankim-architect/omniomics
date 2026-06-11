#!/usr/bin/env python3
"""Does RNA-vs-enhancer-methylation DISAGREEMENT add value on different SUBTYPE AXES?
Tests three real breast-cancer axes with axis-appropriate pole priors, enhancer methylation for
all samples. Repeated CV + paired Wilcoxon for (a) multi-omics gain and (b) disagreement gain."""
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
M=pd.read_csv(os.path.join(MT,"enh_meth_all.tsv.gz"),sep="\t",index_col=0).astype("float32")
meth=me.aggregate(M, os.path.join(MT,"pg_enhancer.tsv")); del M
print(f"[axes] enhancer methylation: {meth.shape[0]} genes x {meth.shape[1]} samples")

MARK={"ER":["ESR1","GATA3","FOXA1","XBP1","PGR"],
      "PROLIF":["MKI67","CCNB1","CCNE1","BUB1","AURKA","MYBL2","CDK1","E2F1","FOXM1"],
      "HER2":["ERBB2","GRB7","STARD3","PGAP3","MIEN1"],
      "BASAL":["KRT5","KRT14","KRT17","FOXC1","MIA","SFRP1","EGFR"]}
AXES={
 "HER2-vs-Luminal":  ("cohort_v4.tsv","HER2",   {"HER2":["HER2"],"Lum":["ER"]}),
 "LumA-vs-LumB":     ("cohort_v2.tsv","LumB",   {"ER":["ER"],"Prolif":["PROLIF"]}),
 "Basal-vs-Luminal": ("cohort_v3.tsv","Basal",  {"Basal":["BASAL"],"Lum":["ER"]}),
}
gmt={k:v for k,v in MARK.items()}

def rep_cv(Xdict,y,rep=10):
    out={k:[] for k in Xdict}
    for r in range(rep):
        sp=list(StratifiedKFold(5,shuffle=True,random_state=r).split(np.zeros(len(y)),y))
        for k,X in Xdict.items():
            a=[]
            for tr,te in sp:
                sc=StandardScaler().fit(X[tr]); clf=LogisticRegression(max_iter=3000,C=1.0).fit(sc.transform(X[tr]),y[tr])
                a.append(roc_auc_score(y[te],clf.predict_proba(sc.transform(X[te]))[:,1]))
            out[k].append(np.mean(a))
    return {k:np.array(v) for k,v in out.items()}

rows=[]
for axis,(cf,pos,poles) in AXES.items():
    coh=pd.read_csv(os.path.join(BR,cf),sep="\t").set_index("sample_id")
    lab=coh[coh["group"].isin([pos]+[g for g in coh["group"].unique() if g!=pos])]["group"]
    # binary: pos vs the other dominant class
    classes=lab.value_counts().index.tolist(); other=[c for c in classes if c!=pos][0]
    lab=lab[lab.isin([pos,other])]
    S=[s for s in lab.index if s in rna.columns and s in meth.columns]
    y=(lab.reindex(S)==pos).astype(int).values
    G=sorted(set(rna.index)&set(meth.index)); R=rna.loc[G,S]; Me=meth.loc[G,S]
    D=mo.dmoi_representation(R,Me,poles,gmt)
    pk=list(poles)
    feats={"RNA pole":D[[f"rna_{p}" for p in pk]].values,
           "RNA+meth":D[[f"{a}_{p}" for a in ["rna","meth"] for p in pk]].values,
           "RNA+meth+disagree":D.values}
    Rr=rep_cv(feats,y)
    mo_gain=Rr["RNA+meth"].mean()-Rr["RNA pole"].mean(); p_mo=stats.wilcoxon(Rr["RNA+meth"],Rr["RNA pole"]).pvalue
    dis_gain=Rr["RNA+meth+disagree"].mean()-Rr["RNA+meth"].mean(); p_dis=stats.wilcoxon(Rr["RNA+meth+disagree"],Rr["RNA+meth"]).pvalue
    rows.append({"axis":axis,"n":len(S),"RNA_pole":Rr["RNA pole"].mean(),"RNA+meth":Rr["RNA+meth"].mean(),
                 "+disagree":Rr["RNA+meth+disagree"].mean(),"multiomics_gain":mo_gain,"mo_p":p_mo,
                 "disagree_gain":dis_gain,"disagree_p":p_dis})
res=pd.DataFrame(rows).set_index("axis"); res.to_csv(os.path.join(ROOT,"disagreement_axes.csv"))
pd.set_option("display.width",160)
print("\n[verdict] disagreement value across subtype axes (10x5 CV):")
print(res.round(4).to_string())

fig,ax=plt.subplots(figsize=(9,4.3)); x=np.arange(len(res)); w=0.27
ax.bar(x-w,res["RNA_pole"],w,label="RNA pole",color="#999")
ax.bar(x,  res["RNA+meth"],w,label="RNA+enhMeth",color="#bb8844")
ax.bar(x+w,res["+disagree"],w,label="+disagreement",color="#2E5A87")
ax.set_xticks(x); ax.set_xticklabels(res.index,fontsize=9); ax.set_ylim(0.5,1.02)
ax.set_ylabel("AUROC (10x5 CV)"); ax.set_title("Multi-omics & disagreement value by subtype axis")
for i,(b,d,pm,pd_) in enumerate(zip(res["multiomics_gain"],res["disagree_gain"],res["mo_p"],res["disagree_p"])):
    ax.text(i, res["+disagree"].iloc[i]+0.01, f"mo+{b:.3f}{'*' if pm<0.05 else ''}\ndis+{d:.3f}{'*' if pd_<0.05 else ''}", ha="center",fontsize=7)
ax.legend(fontsize=8); plt.tight_layout(); plt.savefig(os.path.join(FIG,"disagreement_axes.png"),dpi=140); plt.close()
print("[axes] saved figure + metrics  (* = paired p<0.05)")
