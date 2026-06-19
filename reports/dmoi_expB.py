import numpy as np, pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from scipy.stats import spearmanr, wilcoxon

df=pd.read_csv("/sessions/sleepy-blissful-allen/mnt/molcell_epigenetics/tri_omic_genes.csv")
df=df.dropna(subset=["dMeth_WWDminusWT","RNA_log2FC_WWDvsWT"]).copy()
def zz(s): s=s.astype(float); return ((s-s.mean())/s.std()).values
zO=zz(df.dOcc_WWDminusWT); zM=zz(df.dMeth_WWDminusWT)
y=df.RNA_log2FC_WWDvsWT.values.astype(float)
main=np.column_stack([zO,zM]); diff=np.column_stack([zO,zM,zO-zM]); prod=np.column_stack([zO,zM,zO*zM])
def oof(makem,X,rep=10):
    r2=[];rho=[]
    for r in range(rep):
        pred=np.zeros(len(y))
        for tr,te in KFold(5,shuffle=True,random_state=r).split(X):
            pred[te]=makem().fit(X[tr],y[tr]).predict(X[te])
        r2.append(r2_score(y,pred)); rho.append(spearmanr(pred,y).correlation)
    return np.array(r2),np.array(rho)
lin=lambda:make_pipeline(StandardScaler(),Ridge(alpha=1.0))
gbt=lambda:GradientBoostingRegressor(n_estimators=120,max_depth=3,random_state=0)
res={"meth only (1, linear)":oof(lin,zM.reshape(-1,1)),
     "occ+meth (2, linear)":oof(lin,main),
     "+difference (3, linear)":oof(lin,diff),
     "+product (3, linear)":oof(lin,prod),
     "GBT occ+meth (nonlinear)":oof(gbt,main)}
base=res["occ+meth (2, linear)"][0]
print(f"GSE57577 gene-level: predict WWD-vs-WT log2FC from occupancy/methylation deltas  (n={len(y)}, 10x5 OOF CV)")
for k in ["meth only (1, linear)","occ+meth (2, linear)","+difference (3, linear)","+product (3, linear)","GBT occ+meth (nonlinear)"]:
    r2,rho=res[k]; p=1.0 if k=="occ+meth (2, linear)" else wilcoxon(res[k][0],base).pvalue
    tag="" if k=="occ+meth (2, linear)" else ("  <-- GAIN" if (r2.mean()-base.mean()>0 and p<0.05) else "  (ns/neg)")
    print(f"  {k:28s} R2={r2.mean():+.4f}  Spearman(pred,y)={rho.mean():+.3f}  dR2={r2.mean()-base.mean():+.4f}  p={p:.1e}{tag}")
