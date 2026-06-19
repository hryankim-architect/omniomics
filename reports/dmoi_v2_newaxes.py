import numpy as np, pandas as pd, sys
sys.path.insert(0,"/sessions/sleepy-blissful-allen/mnt/molcell_epigenetics/omniomics-prototype")
from omniomics import multiomics as mo
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedShuffleSplit, StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, brier_score_loss
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
r15=pd.read_csv(f"{OUT}/rna1500.tsv",sep="\t",index_col=0)
S=[s for s in S0 if s in r15.columns]; y=np.array([1 if lab[s]=="LumB" else 0 for s in S])
pg=[g for g in poleg if g in rsub.index and g in meth.index]
rna_gl=rsub.loc[pg]; meth_gl=meth.loc[pg]
Xv2=mo.dmoi_v2_representation(rna_gl[S], meth_gl[S], POLES, GMT).loc[S].values
Xr=r15[S].T.values
def Vfeat(kept):
    k=[g for g in pg if g in kept]
    return mo.dmoi_v2_representation(rna_gl.loc[k][S], meth_gl.loc[k][S], POLES, GMT).loc[S].values

# ---- B: calibration (Brier, lower=better), OOF 5-fold ----
cv=StratifiedKFold(5,shuffle=True,random_state=0)
pR=cross_val_predict(make_pipeline(StandardScaler(),LogisticRegression(max_iter=4000,C=0.5)),Xr,y,cv=cv,method="predict_proba")[:,1]
pV=cross_val_predict(GradientBoostingClassifier(n_estimators=80,max_depth=2,random_state=0),Xv2,y,cv=cv,method="predict_proba")[:,1]
print(f"[B calibration] Brier  RNA1500={brier_score_loss(y,pR):.4f}   DMOI v2={brier_score_loss(y,pV):.4f}  (lower=better)")
print(f"[B calibration] AUROC(OOF) RNA1500={roc_auc_score(y,pR):.4f}  v2={roc_auc_score(y,pV):.4f}")

# ---- A: test-time gene dropout robustness ----
(tr,te),=StratifiedShuffleSplit(1,test_size=0.30,random_state=0).split(np.zeros(len(y)),y)
scR=StandardScaler().fit(Xr[tr]); mR=LogisticRegression(max_iter=4000,C=0.5).fit(scR.transform(Xr[tr]),y[tr])
mV=GradientBoostingClassifier(n_estimators=80,max_depth=2,random_state=0).fit(Xv2[tr],y[tr])
rng=np.random.default_rng(0); trmean=Xr[tr].mean(0)
print("[A gene-dropout]  f     RNA1500  DMOIv2")
for f in [0.0,0.2,0.4,0.6,0.8]:
    aR=[];aV=[]
    for _ in range(6):
        keep=rng.random(Xr.shape[1])>=f
        Xd=Xr[te].copy(); Xd[:,~keep]=trmean[~keep]
        aR.append(roc_auc_score(y[te],mR.predict_proba(scR.transform(Xd))[:,1]))
        keptg=set(g for g in pg if rng.random()>=f)
        if len(keptg)<1: keptg=set(list(pg)[:1])
        Xvd=Vfeat(keptg)[te]
        aV.append(roc_auc_score(y[te],mV.predict_proba(Xvd)[:,1]))
    print(f"   {f:.1f}   {np.mean(aR):.4f}   {np.mean(aV):.4f}")

# ---- C: complementarity on RNA-ambiguous test cases ----
probR_te=mR.predict_proba(scR.transform(Xr[te]))[:,1]
amb=(probR_te>0.35)&(probR_te<0.65)
if amb.sum()>=15:
    yt=y[te][amb]
    aR=roc_auc_score(yt,probR_te[amb])
    aV=roc_auc_score(yt,mV.predict_proba(Xv2[te][amb])[:,1])
    print(f"[C ambiguous cases] n={int(amb.sum())} where RNA prob in [0.35,0.65]:  RNA AUROC={aR:.3f}  v2 AUROC={aV:.3f}")
else:
    print(f"[C ambiguous cases] only {int(amb.sum())} ambiguous (too few to judge)")
