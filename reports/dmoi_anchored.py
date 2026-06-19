import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit, KFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score
from scipy.stats import wilcoxon
D="/sessions/sleepy-blissful-allen/mnt/AI/dmoi-brca-poc/data/tcga_brca"; OUT="/sessions/sleepy-blissful-allen/mnt/outputs"
coh=pd.read_csv(f"{D}/cohort_v2.tsv",sep="\t"); coh=coh[(coh.group.isin(["LumA","LumB"]))&(coh.has_rna)&(coh.has_meth)]
lab=dict(zip(coh.sample_id,coh.group))
r15=pd.read_csv(f"{OUT}/rna1500.tsv",sep="\t",index_col=0); msub=pd.read_csv(f"{OUT}/meth_sub.tsv",sep="\t",index_col=0)
S=[s for s in r15.columns if s in lab and s in msub.columns]; y=np.array([1 if lab[s]=="LumB" else 0 for s in S])
Xr=r15[S].T.values; Mm=msub[S].T.astype(float); Mm=Mm.fillna(Mm.mean()); Xm=Mm.values
def logit(p): p=np.clip(p,1e-4,1-1e-4); return np.log(p/(1-p))
rna_lr=lambda: make_pipeline(StandardScaler(),LogisticRegression(max_iter=5000,C=0.5))
meth_ridge=lambda: make_pipeline(StandardScaler(),Ridge(alpha=20.0))
betas=[0,0.25,0.5,1,2,4]
aAnchor=[];aGated=[];chosen=[]
for r in range(10):
    (tr,te),=StratifiedShuffleSplit(1,test_size=0.30,random_state=r).split(np.zeros(len(y)),y); ytr=y[tr]
    pR_in=cross_val_predict(rna_lr(),Xr[tr],ytr,cv=StratifiedKFold(5,shuffle=True,random_state=r),method="predict_proba")[:,1]
    resid=ytr-pR_in
    mc_in=cross_val_predict(meth_ridge(),Xm[tr],resid,cv=KFold(5,shuffle=True,random_state=r))
    base_in=logit(pR_in)
    b=betas[int(np.argmax([roc_auc_score(ytr, base_in+bb*mc_in) for bb in betas]))]; chosen.append(b)
    pR_te=rna_lr().fit(Xr[tr],ytr).predict_proba(Xr[te])[:,1]
    mc_te=meth_ridge().fit(Xm[tr],resid).predict(Xm[te])
    aAnchor.append(roc_auc_score(y[te],pR_te)); aGated.append(roc_auc_score(y[te], logit(pR_te)+b*mc_te))
aAnchor=np.array(aAnchor);aGated=np.array(aGated)
print(f"n={len(S)}  10x (70/30) outer; beta chosen on train-inner-OOF (beta=0 allowed -> defaults to anchor)")
print(f"  RNA anchor alone           AUROC={aAnchor.mean():.4f}")
print(f"  RNA-anchored GATED + meth   AUROC={aGated.mean():.4f}   dVsAnchor={(aGated-aAnchor).mean():+.4f}  p={wilcoxon(aGated,aAnchor).pvalue:.2e}")
print(f"  chosen beta per repeat: {chosen}")
print(f"  gated >= anchor in {int((aGated>=aAnchor-1e-9).sum())}/10 repeats")
