import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.model_selection import StratifiedShuffleSplit, StratifiedKFold, KFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score
from scipy.stats import wilcoxon
D="/sessions/sleepy-blissful-allen/mnt/AI/dmoi-brca-poc/data/tcga_brca"; OUT="/sessions/sleepy-blissful-allen/mnt/outputs"
cl=pd.read_csv(f"{D}/BRCA_clinicalMatrix.tsv",sep="\t",usecols=["sampleID","age_at_initial_pathologic_diagnosis"])
cl["age"]=pd.to_numeric(cl["age_at_initial_pathologic_diagnosis"],errors="coerce"); cl=cl.dropna(subset=["age"])
age=dict(zip(cl.sampleID,cl.age))
r15=pd.read_csv(f"{OUT}/rna1500.tsv",sep="\t",index_col=0)
mgw=pd.read_csv(f"{OUT}/meth_gw.tsv",sep="\t",index_col=0)
S=[s for s in r15.columns if s in mgw.columns and s in age]
ages=np.array([age[s] for s in S]); med=np.median(ages); y=(ages>med).astype(int)
Xr=r15[S].T.values
Mg=mgw[S].T.astype(float); Mg=Mg.fillna(Mg.mean()).fillna(0.0); Xm=Mg.values
def logit(p): p=np.clip(p,1e-4,1-1e-4); return np.log(p/(1-p))
rna_lr=lambda: make_pipeline(StandardScaler(),LogisticRegression(max_iter=5000,C=0.3))
meth_lr=lambda: make_pipeline(StandardScaler(),LogisticRegression(max_iter=5000,C=0.3))
meth_ridge=lambda: make_pipeline(StandardScaler(),Ridge(alpha=50.0))
betas=[0,0.25,0.5,1,2,4,8]
aR=[];aM=[];aG=[];bs=[]
for r in range(8):
    (tr,te),=StratifiedShuffleSplit(1,test_size=0.30,random_state=r).split(np.zeros(len(y)),y); ytr=y[tr]
    pRin=cross_val_predict(rna_lr(),Xr[tr],ytr,cv=StratifiedKFold(5,shuffle=True,random_state=r),method="predict_proba")[:,1]
    resid=ytr-pRin
    mcin=cross_val_predict(meth_ridge(),Xm[tr],resid,cv=KFold(5,shuffle=True,random_state=r))
    base=logit(pRin)
    b=betas[int(np.argmax([roc_auc_score(ytr,base+bb*mcin) for bb in betas]))]; bs.append(b)
    pRte=rna_lr().fit(Xr[tr],ytr).predict_proba(Xr[te])[:,1]
    mcte=meth_ridge().fit(Xm[tr],resid).predict(Xm[te])
    pMte=meth_lr().fit(Xm[tr],ytr).predict_proba(Xm[te])[:,1]
    aR.append(roc_auc_score(y[te],pRte)); aM.append(roc_auc_score(y[te],pMte)); aG.append(roc_auc_score(y[te],logit(pRte)+b*mcte))
aR=np.array(aR);aM=np.array(aM);aG=np.array(aG)
print(f"AGE endpoint (older-than-median)  n={len(S)}  (pos={int(y.sum())})  8x(70/30)")
print(f"  RNA alone (anchor)             AUROC={aR.mean():.4f}")
print(f"  Methylation alone (3036 CpG)   AUROC={aM.mean():.4f}   <-- epigenetic-clock signal")
print(f"  RNA-anchored, GATED + meth     AUROC={aG.mean():.4f}   dVsRNA={(aG-aR).mean():+.4f}  p={wilcoxon(aG,aR).pvalue:.2e}")
print(f"  chosen beta per repeat: {bs}   (beta>0 => methylation earns its place)")
