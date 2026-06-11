#!/usr/bin/env python3
"""Where does the RNA-vs-methylation DISAGREEMENT signal earn value?
On LumA-vs-LumB with enhancer methylation, degrade RNA information (inject noise into the RNA
pole scores) and ask whether enhancer-methylation + disagreement compensates. Hypothesis: the
multi-omics/disagreement value grows as RNA becomes information-poor."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from omniomics import loaders, methylation as me, multiomics as mo
from omniomics import config
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BR=config.brca_tcga_dir()
MT=config.cache_dir(); GMT=config.hallmark_gmt()
ROOT=os.path.dirname(__file__); FIG=os.path.join(ROOT,"figures"); rng=np.random.default_rng(0)

coh=pd.read_csv(os.path.join(BR,"cohort_v2.tsv"),sep="\t").set_index("sample_id")
lab=coh[coh["group"].isin(["LumA","LumB"])]["group"]
rna=loaders.gene_symbol_index(loaders.load_matrix(os.path.join(BR,"HiSeqV2.gz"),gene_col=0))
M=pd.read_csv(os.path.join(MT,"ctx_meth_lumab.tsv.gz"),sep="\t",index_col=0).astype("float32")
meth_enh=me.aggregate(M, os.path.join(MT,"pg_enhancer.tsv")); del M
S=[s for s in lab.index if s in rna.columns and s in meth_enh.columns]
y=(lab.reindex(S)=="LumB").astype(int).values
G=sorted(set(rna.index)&set(meth_enh.index)); rna=rna.loc[G,S]; meth=meth_enh.loc[G,S]
gmt=mo.load_gmt(GMT)
POLE={"LumA":["HALLMARK_ESTROGEN_RESPONSE_EARLY","HALLMARK_ESTROGEN_RESPONSE_LATE"],
      "LumB":["HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MYC_TARGETS_V1"]}
D=mo.dmoi_representation(rna, meth, POLE, gmt)         # rna_/meth_/disagree_ per pole
rna_pole=D[[f"rna_{p}" for p in POLE]].values
meth_pole=D[[f"meth_{p}" for p in POLE]].values
print(f"[regime] {len(S)} samples; degrading RNA, testing methylation/disagreement compensation")

def rcv(X,y,rep=10):
    a=[]
    for r in range(rep):
        for tr,te in StratifiedKFold(5,shuffle=True,random_state=r).split(np.zeros(len(y)),y):
            sc=StandardScaler().fit(X[tr]); clf=LogisticRegression(max_iter=3000,C=1.0).fit(sc.transform(X[tr]),y[tr])
            a.append(roc_auc_score(y[te],clf.predict_proba(sc.transform(X[te]))[:,1]))
    return np.mean(a)

sigmas=[0.0,0.5,1.0,2.0,3.0]; rows=[]
zr=StandardScaler().fit_transform(rna_pole)
for s in sigmas:
    rdeg = zr + rng.normal(0,s,size=zr.shape)                  # noisy RNA
    disagree = rdeg - StandardScaler().fit_transform(meth_pole)
    a_rna  = rcv(rdeg, y)
    a_rm   = rcv(np.hstack([rdeg, meth_pole]), y)
    a_dmoi = rcv(np.hstack([rdeg, meth_pole, disagree]), y)
    rows.append({"RNA_noise_sigma":s,"RNA_only":a_rna,"RNA+meth":a_rm,"RNA+meth+disagree":a_dmoi,
                 "multiomics_gain":a_dmoi-a_rna})
res=pd.DataFrame(rows).set_index("RNA_noise_sigma"); res.to_csv(os.path.join(ROOT,"disagreement_regime.csv"))
print("\n[verdict] AUROC vs RNA degradation:")
print(res.round(3).to_string())

fig,ax=plt.subplots(1,2,figsize=(11,4.3))
ax[0].plot(res.index,res["RNA_only"],"o-",label="RNA only")
ax[0].plot(res.index,res["RNA+meth"],"s-",label="RNA + enhancer-meth")
ax[0].plot(res.index,res["RNA+meth+disagree"],"^-",label="+ disagreement (DMOI)")
ax[0].set_xlabel("RNA noise σ (information loss →)"); ax[0].set_ylabel("LumA-vs-LumB AUROC"); ax[0].legend(fontsize=8); ax[0].set_title("Methylation compensates as RNA degrades")
ax[1].plot(res.index,res["multiomics_gain"],"D-",color="#2E5A87"); ax[1].axhline(0,ls="--",c="k",lw=.6)
ax[1].set_xlabel("RNA noise σ"); ax[1].set_ylabel("multi-omics gain (DMOI − RNA)"); ax[1].set_title("Disagreement/methylation value grows when RNA is poor")
plt.tight_layout(); plt.savefig(os.path.join(FIG,"disagreement_regime.png"),dpi=140); plt.close()
print("[regime] saved figure + metrics")
