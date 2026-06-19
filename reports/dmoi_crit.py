import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from scipy import stats
D="/sessions/sleepy-blissful-allen/mnt/AI/dmoi-brca-poc/data/tcga_brca"
POLES={"ER":["ESR1","GATA3","FOXA1","XBP1","PGR"],
       "PROLIF":["MKI67","CCNB1","CCNE1","BUB1","AURKA","MYBL2","CDK1","E2F1","FOXM1"]}
poleg=[g for v in POLES.values() for g in v]
coh=pd.read_csv(f"{D}/cohort_v2.tsv",sep="\t")
coh=coh[(coh.group.isin(["LumA","LumB"]))&(coh.has_rna)&(coh.has_meth)]
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
methgene={g: meth.loc[[p for p in g2p[g] if p in meth.index]].mean(axis=0) for g in poleg if any(p in meth.index for p in g2p[g])}
METH=pd.DataFrame(methgene)
S=[s for s in rna.columns if s in lab and s in METH.index]
y=np.array([1 if lab[s]=="LumB" else 0 for s in S])
RNA=rna[S].T
METH=METH.loc[S]
def zc(df):
    df=df.astype(float); df=df.fillna(df.mean()); return (df-df.mean())/(df.std(ddof=0)+1e-9)
zR=zc(RNA); zM=zc(-METH)
def pole(z,genes):
    c=[g for g in genes if g in z.columns]; return z[c].mean(axis=1)
F={}
for p,genes in POLES.items(): F[f"rna_{p}"]=pole(zR,genes); F[f"meth_{p}"]=pole(zM,genes)
F=pd.DataFrame(F); pk=list(POLES)
RNAm=F[[f"rna_{p}" for p in pk]+[f"meth_{p}" for p in pk]].values
dis =np.column_stack([F[f"rna_{p}"]-F[f"meth_{p}"] for p in pk])
prod=np.column_stack([F[f"rna_{p}"]*F[f"meth_{p}"] for p in pk])
feats={"RNA+meth (4, linear)":RNAm,"+disagree=difference (6)":np.column_stack([RNAm,dis]),"+concordance=product (6)":np.column_stack([RNAm,prod])}
def cv(X,gbt=False,rep=12):
    o=[]
    for r in range(rep):
        a=[]
        for tr,te in StratifiedKFold(5,shuffle=True,random_state=r).split(X,y):
            if gbt:
                clf=GradientBoostingClassifier(n_estimators=80,max_depth=2,random_state=r).fit(X[tr],y[tr]); pr=clf.predict_proba(X[te])[:,1]
            else:
                sc=StandardScaler().fit(X[tr]); clf=LogisticRegression(max_iter=5000).fit(sc.transform(X[tr]),y[tr]); pr=clf.predict_proba(sc.transform(X[te]))[:,1]
            a.append(roc_auc_score(y[te],pr))
        o.append(np.mean(a))
    return np.array(o)
R={k:cv(X) for k,X in feats.items()}
R["GBT on RNA+meth (nonlinear)"]=cv(RNAm,gbt=True)
base=R["RNA+meth (4, linear)"]
print(f"n={len(S)}  LumA={int((1-y).sum())}  LumB={int(y.sum())}   (12x5 CV)")
for k in ["RNA+meth (4, linear)","+disagree=difference (6)","+concordance=product (6)","GBT on RNA+meth (nonlinear)"]:
    p=1.0 if k=="RNA+meth (4, linear)" else stats.wilcoxon(R[k],base).pvalue
    d=(R[k]-base).mean()
    tag="" if k=="RNA+meth (4, linear)" else ("  <-- GAIN" if (d>0 and p<0.05) else "  (ns/neg)")
    print(f"  {k:34s} AUROC={R[k].mean():.4f}  d={d:+.4f}  p={p:.1e}{tag}")
