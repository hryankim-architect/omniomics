import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
D="/sessions/sleepy-blissful-allen/mnt/AI/dmoi-brca-poc/data/tcga_brca"
POLES={"ER":["ESR1","GATA3","FOXA1","XBP1","PGR"],"PROLIF":["MKI67","CCNB1","CCNE1","BUB1","AURKA","MYBL2","CDK1","E2F1","FOXM1"]}
poleg=[g for v in POLES.values() for g in v]
coh=pd.read_csv(f"{D}/cohort_v2.tsv",sep="\t"); coh=coh[(coh.group.isin(["LumA","LumB"]))&(coh.has_rna)&(coh.has_meth)]
lab=dict(zip(coh.sample_id,coh.group))
rna=pd.read_csv("rna_sub.tsv",sep="\t",index_col=0); meth=pd.read_csv("meth_sub.tsv",sep="\t",index_col=0)
g2p={g:[] for g in poleg}
with open(f"{D}/hm450_probemap.tsv") as fh:
    next(fh)
    for line in fh:
        f=line.split("\t")
        if len(f)>=2 and f[1]!=".":
            for g in f[1].split(","):
                if g in g2p: g2p[g].append(f[0])
S=[s for s in rna.columns if s in lab and s in meth.columns]
y=np.array([1 if lab[s]=="LumB" else 0 for s in S])
def impute(arr):  # fill NaN with per-probe(row) mean
    arr=arr.copy(); rm=np.nanmean(arr,axis=1,keepdims=True); ii=np.where(np.isnan(arr)); arr[ii]=np.take(rm,ii[0]); return np.nan_to_num(arr,nan=0.0)
probeM={g: impute(meth.loc[[p for p in g2p[g] if p in meth.index],S].astype(float).values) for g in poleg if any(p in meth.index for p in g2p[g])}
RNA=rna[S].T.astype(float); zR=(RNA-RNA.mean())/(RNA.std(ddof=0)+1e-9)
gstd=float(np.std(np.concatenate([m.ravel() for m in probeM.values()])))
def pole_rna(genes): c=[g for g in genes if g in zR.columns]; return zR[c].mean(axis=1).values
rnaP=np.column_stack([pole_rna(POLES[p]) for p in POLES])
rng=np.random.default_rng(0)
def features(alpha, weighted):
    methg={}; rel={}
    for g,Mp in probeM.items():
        Mn=Mp+rng.normal(0,alpha*gstd,Mp.shape) if alpha>0 else Mp
        methg[g]=Mn.mean(0)
        rel[g]=(np.nanmean(np.corrcoef(Mn)[np.triu_indices(len(Mn),1)]) if len(Mn)>=2 else np.nan)
    med=np.nanmedian([v for v in rel.values() if v==v]); rel={g:(rel[g] if rel[g]==rel[g] else med) for g in methg}
    METH=pd.DataFrame(methg,index=S); zM=(-METH-(-METH).mean())/((-METH).std(ddof=0)+1e-9)
    def pm(genes):
        c=[g for g in genes if g in zM.columns]
        if not c: return np.zeros(len(S))
        if not weighted: return zM[c].mean(1).values
        w=np.array([max(rel[g],0.0)+1e-3 for g in c]); return (zM[c].values*w).sum(1)/w.sum()
    methP=np.column_stack([pm(POLES[p]) for p in POLES])
    intP =np.column_stack([rnaP[:,i]*methP[:,i] for i in range(len(POLES))])
    return np.column_stack([rnaP,methP,intP])
def cv(X,rep=6):
    return np.mean([np.mean([roc_auc_score(y[te],GradientBoostingClassifier(n_estimators=80,max_depth=2,random_state=r).fit(X[tr],y[tr]).predict_proba(X[te])[:,1]) for tr,te in StratifiedKFold(5,shuffle=True,random_state=r).split(X,y)]) for r in range(rep)])
rf=np.mean([np.mean([roc_auc_score(y[te],LogisticRegression(max_iter=4000).fit(StandardScaler().fit_transform(rnaP[tr]),y[tr]).predict_proba(StandardScaler().fit(rnaP[tr]).transform(rnaP[te]))[:,1]) for tr,te in StratifiedKFold(5,shuffle=True,random_state=r).split(rnaP,y)]) for r in range(6)])
alphas=[0,1,2,4,8]
eq=[cv(features(a,False)) for a in alphas]; wt=[cv(features(a,True)) for a in alphas]
print("alpha (noise sd / meth sd):", alphas)
print("equal-weight fusion  :", [f"{v:.4f}" for v in eq])
print("reliability-weighted :", [f"{v:.4f}" for v in wt])
print(f"RNA-only floor       : {rf:.4f}")
fig,ax=plt.subplots(figsize=(6.4,4.2),dpi=200)
ax.plot(alphas,eq,"o-",color="#C0504D",lw=2.2,label="equal-weight fusion")
ax.plot(alphas,wt,"s-",color="#1F3B5B",lw=2.2,label="reliability-weighted (DMOI v2)")
ax.axhline(rf,ls="--",color="#888",lw=1.2,label=f"RNA-only floor ({rf:.3f})")
ax.set_xlabel("methylation-layer corruption  (noise sd / methylation sd)")
ax.set_ylabel("LumA-vs-LumB AUROC (6x5 CV)")
ax.set_title("Reliability-gating protects against a degraded layer")
ax.legend(fontsize=8,frameon=False); ax.spines[["top","right"]].set_visible(False); ax.grid(axis="y",color="#eee")
fig.tight_layout(); fig.savefig("/sessions/sleepy-blissful-allen/mnt/molcell_epigenetics/omniomics-prototype/reports/figs/dmoi_v2_robustness.png",bbox_inches="tight")
print("saved figs/dmoi_v2_robustness.png")
