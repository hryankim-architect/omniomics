import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.model_selection import StratifiedShuffleSplit, StratifiedKFold, KFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score
from scipy.stats import wilcoxon
OUT="/sessions/sleepy-blissful-allen/mnt/outputs"
r15=pd.read_csv(f"{OUT}/rna1500.tsv",sep="\t",index_col=0)
mgw=pd.read_csv(f"{OUT}/meth_gw.tsv",sep="\t",index_col=0)
S=[s for s in r15.columns if s in mgw.columns]
Mg=mgw[S].T.astype(float); Mg=Mg.fillna(Mg.mean()).fillna(0.0)
cps=list(Mg.columns); half=len(cps)//2
A=cps[:half]; B=cps[half:]                      # disjoint CpG sets
score=Mg[A].mean(1).values                      # methylation phenotype from set A (RNA cannot see it)
y=(score>np.median(score)).astype(int)
Xr=r15[S].T.values; Xb=Mg[B].values             # predictors: RNA, and DISJOINT methylation set B
def logit(p): p=np.clip(p,1e-4,1-1e-4); return np.log(p/(1-p))
rl=lambda: make_pipeline(StandardScaler(),LogisticRegression(max_iter=5000,C=0.3))
bl=lambda: make_pipeline(StandardScaler(),LogisticRegression(max_iter=5000,C=0.3))
br=lambda: make_pipeline(StandardScaler(),Ridge(alpha=50.0))
betas=[0,0.25,0.5,1,2,4,8]
aR=[];aB=[];aG=[];bs=[]
for r in range(8):
    (tr,te),=StratifiedShuffleSplit(1,test_size=0.30,random_state=r).split(np.zeros(len(y)),y); ytr=y[tr]
    pin=cross_val_predict(rl(),Xr[tr],ytr,cv=StratifiedKFold(5,shuffle=True,random_state=r),method="predict_proba")[:,1]
    resid=ytr-pin
    mcin=cross_val_predict(br(),Xb[tr],resid,cv=KFold(5,shuffle=True,random_state=r))
    base=logit(pin)
    b=betas[int(np.argmax([roc_auc_score(ytr,base+bb*mcin) for bb in betas]))]; bs.append(b)
    pRte=rl().fit(Xr[tr],ytr).predict_proba(Xr[te])[:,1]
    mcte=br().fit(Xb[tr],resid).predict(Xb[te])
    pBte=bl().fit(Xb[tr],ytr).predict_proba(Xb[te])[:,1]
    aR.append(roc_auc_score(y[te],pRte)); aB.append(roc_auc_score(y[te],pBte)); aG.append(roc_auc_score(y[te],logit(pRte)+b*mcte))
aR=np.array(aR);aB=np.array(aB);aG=np.array(aG)
print(f"POSITIVE CONTROL: methylation-defined endpoint (mean of CpG set A)  n={len(S)} 8x(70/30)")
print(f"  RNA alone (anchor; cannot see the meth axis)  AUROC={aR.mean():.4f}")
print(f"  Methylation set B alone (disjoint from A)     AUROC={aB.mean():.4f}")
print(f"  RNA-anchored, GATED + methB                   AUROC={aG.mean():.4f}   dVsRNA={(aG-aR).mean():+.4f}  p={wilcoxon(aG,aR).pvalue:.2e}")
print(f"  chosen beta per repeat: {bs}   (beta>0 => the gate engages and gains)")
