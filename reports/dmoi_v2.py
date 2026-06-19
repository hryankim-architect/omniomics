import numpy as np, pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from scipy.stats import wilcoxon
D="/sessions/sleepy-blissful-allen/mnt/AI/dmoi-brca-poc/data/tcga_brca"
POLES={"ER":["ESR1","GATA3","FOXA1","XBP1","PGR"],
       "PROLIF":["MKI67","CCNB1","CCNE1","BUB1","AURKA","MYBL2","CDK1","E2F1","FOXM1"]}
poleg=[g for v in POLES.values() for g in v]
coh=pd.read_csv(f"{D}/cohort_v2.tsv",sep="\t"); coh=coh[(coh.group.isin(["LumA","LumB"]))&(coh.has_rna)&(coh.has_meth)]
lab=dict(zip(coh.sample_id,coh.group))
rna=pd.read_csv("rna_sub.tsv",sep="\t",index_col=0)
meth=pd.read_csv("meth_sub.tsv",sep="\t",index_col=0)
g2p={g:[] for g in poleg}
with open(f"{D}/hm450_probemap.tsv") as fh:
    next(fh)
    for line in fh:
        f=line.split("\t")
        if len(f)<2 or f[1]==".": continue
        for g in f[1].split(","):
            if g in g2p: g2p[g].append(f[0])
S=[s for s in rna.columns if s in lab and s in meth.columns]
y=np.array([1 if lab[s]=="LumB" else 0 for s in S])
# per-gene methylation (mean over probes) + reliability (mean inter-probe corr across samples)
methg={}; rel={}
for g in poleg:
    ps=[p for p in g2p[g] if p in meth.index]
    if not ps: continue
    M=meth.loc[ps,S].astype(float)              # probes x samples
    methg[g]=M.mean(axis=0).values
    if len(ps)>=2:
        C=np.corrcoef(M.values)                 # probe x probe
        iu=np.triu_indices(len(ps),1); rel[g]=np.nanmean(C[iu])
    else: rel[g]=np.nan
med=np.nanmedian([v for v in rel.values() if v==v])
rel={g:(rel[g] if rel.get(g)==rel.get(g) else med) for g in methg}
METH=pd.DataFrame(methg, index=S)               # samples x genes
RNA=rna[S].T                                     # samples x genes
def zc(df): df=df.astype(float).fillna(df.astype(float).mean()); return (df-df.mean())/(df.std(ddof=0)+1e-9)
zR=zc(RNA); zM=zc(-METH)                         # meth_sign=-1
def pole_rna(genes): c=[g for g in genes if g in zR.columns]; return zR[c].mean(axis=1).values
def pole_meth(genes,weighted):
    c=[g for g in genes if g in zM.columns]
    if not c: return np.zeros(len(S))
    if not weighted: return zM[c].mean(axis=1).values
    w=np.array([max(rel[g],0.0)+1e-3 for g in c]); W=zM[c].values*w; return W.sum(1)/w.sum()
pk=list(POLES)
rnaP=np.column_stack([pole_rna(POLES[p]) for p in pk])
methU=np.column_stack([pole_meth(POLES[p],False) for p in pk])
methR=np.column_stack([pole_meth(POLES[p],True)  for p in pk])
diffU=np.column_stack([rnaP[:,i]-methU[:,i] for i in range(len(pk))])
intR =np.column_stack([rnaP[:,i]*methR[:,i] for i in range(len(pk))])
intU =np.column_stack([rnaP[:,i]*methU[:,i] for i in range(len(pk))])
FE={
 "RNA poles only (LR)":            (rnaP, False),
 "orig DMOI: RNA+meth+diff (LR)":  (np.column_stack([rnaP,methU,diffU]), False),
 "v2 no-reliability: +int (GBT)":  (np.column_stack([rnaP,methU,intU]), True),
 "v2 no-interaction: relMeth (GBT)":(np.column_stack([rnaP,methR]), True),
 "DMOI v2: relMeth+int (GBT)":     (np.column_stack([rnaP,methR,intR]), True),
}
def cv(X,gbt,rep=12):
    o=[]
    for r in range(rep):
        a=[]
        for tr,te in StratifiedKFold(5,shuffle=True,random_state=r).split(X,y):
            if gbt: clf=GradientBoostingClassifier(n_estimators=100,max_depth=2,random_state=r).fit(X[tr],y[tr]); pr=clf.predict_proba(X[te])[:,1]
            else: sc=StandardScaler().fit(X[tr]); clf=LogisticRegression(max_iter=5000).fit(sc.transform(X[tr]),y[tr]); pr=clf.predict_proba(sc.transform(X[te]))[:,1]
            a.append(roc_auc_score(y[te],pr))
        o.append(np.mean(a))
    return np.array(o)
R={k:cv(X,g) for k,(X,g) in FE.items()}
base=R["orig DMOI: RNA+meth+diff (LR)"]
print(f"BRCA LumA/B  n={len(S)} (LumB={int(y.sum())})  12x5 CV   [reliability(rel) per gene: "+", ".join(f"{g}:{rel[g]:.2f}" for g in list(rel)[:6])+" ...]")
print(f"  reference: plain ~1500-gene RNA (repo) = 0.940")
for k in FE:
    p=1.0 if k=="orig DMOI: RNA+meth+diff (LR)" else wilcoxon(R[k],base).pvalue
    d=R[k].mean()-base.mean()
    tag="" if k=="orig DMOI: RNA+meth+diff (LR)" else ("  <-- gain vs orig" if (d>0 and p<0.05) else "")
    print(f"  {k:34s} AUROC={R[k].mean():.4f}  dVsOrig={d:+.4f}  p={p:.1e}{tag}")
