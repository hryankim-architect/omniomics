#!/usr/bin/env python3
"""Phase 3c hypothesis test: does LumA-vs-LumB methylation signal live in
ENHANCER / CGI-SHORE probes rather than promoters?
HM450 probes classified by context (promoter TSS+/-1500; CGI-shore 0-2kb from a CpG island;
distal/enhancer >5kb from TSS and >4kb from any CGI). Per-context gene methylation tested."""
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

coh=pd.read_csv(os.path.join(BR,"cohort_v2.tsv"),sep="\t").set_index("sample_id")
lab=coh[coh["group"].isin(["LumA","LumB"])]["group"]
rna=loaders.gene_symbol_index(loaders.load_matrix(os.path.join(BR,"HiSeqV2.gz"),gene_col=0))

print("[ctx] loading context methylation matrix once ...")
M=pd.read_csv(os.path.join(MT,"ctx_meth_lumab.tsv.gz"),sep="\t",index_col=0).astype("float32")
ctx_pg={"promoter":"pg_promoter.tsv","CGI-shore":"pg_shore.tsv","distal/enhancer":"pg_enhancer.tsv"}
meth_ctx={c: me.aggregate(M, os.path.join(MT,pg)) for c,pg in ctx_pg.items()}
del M
for c in meth_ctx: print(f"[ctx] {c}: {meth_ctx[c].shape[0]} genes")

gmt=mo.load_gmt(GMT)
POLE={"LumA":["HALLMARK_ESTROGEN_RESPONSE_EARLY","HALLMARK_ESTROGEN_RESPONSE_LATE"],
      "LumB":["HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MYC_TARGETS_V1"]}
cv=StratifiedKFold(5,shuffle=True,random_state=0)
def auc(X,y): return cross_val_score(LogisticRegression(max_iter=4000,C=1.0),X,y,cv=cv,scoring="roc_auc").mean()

rows=[]
for c, meth in meth_ctx.items():
    S=[s for s in lab.index if s in rna.columns and s in meth.columns]
    y=(lab.reindex(S)=="LumB").astype(int).values
    G=sorted(set(rna.index)&set(meth.index)); R=rna.loc[G,S]; Me=meth.loc[G,S]
    # pole gene sets present in this context
    def pole_genes(p):
        g=set()
        for s in POLE[p]: g|=set(gmt.get(s,[]))
        return [x for x in g if x in Me.index]
    zM=pd.DataFrame(StandardScaler().fit_transform((-Me).T.fillna((-Me).T.mean()).values),index=S,columns=Me.index)
    zR=pd.DataFrame(StandardScaler().fit_transform(R.T.fillna(R.T.mean()).values),index=S,columns=R.index)
    meth_pole=np.column_stack([zM[pole_genes(p)].mean(axis=1).values for p in ["LumA","LumB"]])
    rna_pole =np.column_stack([zR[pole_genes(p)].mean(axis=1).values for p in ["LumA","LumB"]])
    # methylation-alone (top variable context genes)
    var=Me.var(axis=1).sort_values(ascending=False); top=Me.loc[var.index[:1000],S]
    Xm=StandardScaler().fit_transform(top.T.fillna(top.T.mean()).values)
    a_meth_top = auc(Xm,y)
    a_meth_pole= auc(StandardScaler().fit_transform(meth_pole),y)
    a_rna_pole = auc(StandardScaler().fit_transform(rna_pole),y)
    a_combo    = auc(StandardScaler().fit_transform(np.hstack([rna_pole,meth_pole])),y)
    rows.append({"context":c,"n_genes":len(G),"meth_top1000_AUROC":a_meth_top,
                 "meth_pole_AUROC":a_meth_pole,"RNA_pole_AUROC":a_rna_pole,
                 "RNA+meth_pole_AUROC":a_combo,"gain_vs_RNApole":a_combo-a_rna_pole})
res=pd.DataFrame(rows).set_index("context")
res.to_csv(os.path.join(ROOT,"meth_context_auroc.csv"))
print("\n[verdict] LumA-vs-LumB by methylation context:")
print(res.round(3).to_string())

# figure
fig,ax=plt.subplots(figsize=(8,4.2))
x=np.arange(len(res)); w=0.35
ax.bar(x-w/2, res["meth_top1000_AUROC"], w, label="methylation alone (top1000)", color="#bb8844")
ax.bar(x+w/2, res["RNA+meth_pole_AUROC"], w, label="RNA pole + meth pole", color="#2E5A87")
ax.axhline(res["RNA_pole_AUROC"].iloc[0], ls="--", c="k", lw=.8, label="RNA pole only")
ax.axhline(0.5, ls=":", c="gray", lw=.7)
ax.set_xticks(x); ax.set_xticklabels(res.index); ax.set_ylim(0.4,1.0)
ax.set_ylabel("LumA-vs-LumB AUROC (5-fold CV)")
ax.set_title("Does methylation carry subtype signal by genomic context?"); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(FIG,"meth_context_auroc.png"),dpi=140); plt.close()
print("[ctx] saved figure + metrics")
