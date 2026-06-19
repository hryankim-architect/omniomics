import numpy as np, pandas as pd, matplotlib, sys
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0,"/sessions/sleepy-blissful-allen/mnt/molcell_epigenetics/omniomics-prototype")
from omniomics import multiomics as mo
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit
from sklearn.metrics import roc_auc_score
D="/sessions/sleepy-blissful-allen/mnt/AI/dmoi-brca-poc/data/tcga_brca"; OUT="/sessions/sleepy-blissful-allen/mnt/outputs"
POLES={"ER":["ER"],"PROLIF":["PROLIF"]}
GMT={"ER":["ESR1","GATA3","FOXA1","XBP1","PGR"],"PROLIF":["MKI67","CCNB1","CCNE1","BUB1","AURKA","MYBL2","CDK1","E2F1","FOXM1"]}
poleg=[g for v in GMT.values() for g in v]
coh=pd.read_csv(f"{D}/cohort_v2.tsv",sep="\t"); coh=coh[(coh.group.isin(["LumA","LumB"]))&(coh.has_rna)&(coh.has_meth)]
lab=dict(zip(coh.sample_id,coh.group))
rsub=pd.read_csv(f"{OUT}/rna_sub.tsv",sep="\t",index_col=0); msub=pd.read_csv(f"{OUT}/meth_sub.tsv",sep="\t",index_col=0)
g2p={g:[] for g in poleg}
with open(f"{D}/hm450_probemap.tsv") as fh:
    next(fh)
    for line in fh:
        f=line.split("\t")
        if len(f)>=2 and f[1]!=".":
            for g in f[1].split(","):
                if g in g2p: g2p[g].append(f[0])
S0=[s for s in rsub.columns if s in lab and s in msub.columns]
methg={g: msub.loc[[p for p in g2p[g] if p in msub.index],S0].astype(float).mean(0) for g in poleg if any(p in msub.index for p in g2p[g])}
meth=pd.DataFrame(methg).T
D2u=mo.dmoi_v2_representation(rsub[S0], meth[S0], POLES, GMT, reliability=None)
r15=pd.read_csv(f"{OUT}/rna1500.tsv",sep="\t",index_col=0)
S=[s for s in S0 if s in r15.columns]; y=np.array([1 if lab[s]=="LumB" else 0 for s in S])
Xv2=D2u.loc[S].values; Xr=r15[S].T.values
lr=lambda: LogisticRegression(max_iter=4000,C=0.5); gb=lambda: GradientBoostingClassifier(n_estimators=60,max_depth=2,random_state=0)
def cleancv(X,mk,scale,rep=8):
    o=[]
    for r in range(rep):
        a=[]
        for tr,te in StratifiedKFold(5,shuffle=True,random_state=r).split(X,y):
            if scale: sc=StandardScaler().fit(X[tr]); m=mk().fit(sc.transform(X[tr]),y[tr]); p=m.predict_proba(sc.transform(X[te]))[:,1]
            else: m=mk().fit(X[tr],y[tr]); p=m.predict_proba(X[te])[:,1]
            a.append(roc_auc_score(y[te],p))
        o.append(np.mean(a))
    return np.array(o)
accR=cleancv(Xr,lr,True); accV=cleancv(Xv2,gb,False)
sizes=[40,70,120,200,280]; curve={("RNA1500",s):[] for s in sizes}; curve.update({("DMOIv2",s):[] for s in sizes})
for r in range(6):
    (pool,test),=StratifiedShuffleSplit(1,test_size=0.30,random_state=r).split(np.zeros(len(y)),y)
    for s in sizes:
        su=min(s,len(pool)-2)
        (tr0,_),=StratifiedShuffleSplit(1,train_size=su,random_state=100+r).split(np.zeros(len(pool)),y[pool])
        tr=pool[tr0]
        sc=StandardScaler().fit(Xr[tr]); m=lr().fit(sc.transform(Xr[tr]),y[tr]); curve[("RNA1500",s)].append(roc_auc_score(y[test],m.predict_proba(sc.transform(Xr[test]))[:,1]))
        m=gb().fit(Xv2[tr],y[tr]); curve[("DMOIv2",s)].append(roc_auc_score(y[test],m.predict_proba(Xv2[test])[:,1]))
print(f"n={len(S)}  features: RNA1500={Xr.shape[1]}, DMOI v2={Xv2.shape[1]} (named poles)")
print(f"[accuracy 8x5 CV] RNA1500={accR.mean():.4f}±{accR.std():.4f}   DMOI v2={accV.mean():.4f}±{accV.std():.4f}")
print("[data efficiency: AUROC by train size]  n   RNA1500  DMOIv2")
for s in sizes: print(f"  {s:4d}   {np.mean(curve[('RNA1500',s)]):.4f}   {np.mean(curve[('DMOIv2',s)]):.4f}")
fig,ax=plt.subplots(figsize=(6.4,4.3),dpi=200)
for name,col,lbl in [("RNA1500","#C0504D","plain RNA (1500 genes, logistic)"),("DMOIv2","#1F3B5B","DMOI v2 (8 named features, GBT)")]:
    m=np.array([np.mean(curve[(name,s)]) for s in sizes]); sd=np.array([np.std(curve[(name,s)]) for s in sizes])
    ax.plot(sizes,m,"o-",color=col,lw=2.2,label=lbl); ax.fill_between(sizes,m-sd,m+sd,color=col,alpha=0.15)
ax.set_xlabel("training samples"); ax.set_ylabel("LumA-vs-LumB AUROC (held-out 30%)")
ax.set_title("Data efficiency: low-dim interpretable v2 vs high-dim RNA")
ax.legend(fontsize=8,frameon=False,loc="lower right"); ax.spines[["top","right"]].set_visible(False); ax.grid(axis="y",color="#eee")
fig.tight_layout(); fig.savefig("/sessions/sleepy-blissful-allen/mnt/molcell_epigenetics/omniomics-prototype/reports/figs/dmoi_v2_scorecard.png",bbox_inches="tight")
print("saved figs/dmoi_v2_scorecard.png")
