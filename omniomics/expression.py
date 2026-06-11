"""Differential expression: descriptive + paired empirical-Bayes moderated test.
Distilled and verified against the GSE57577 reproduction (scripts 01/04)."""
import numpy as np, pandas as pd
from scipy.special import digamma, polygamma
from scipy import stats

NOISE = [r"^Snor",r"^Snar",r"^Mir\d",r"^Mirlet",r"^n-R5s",r"rRNA$",r"^U\d+$",r"^5S",
         r"^5_8S",r"^7SK",r"^Rny",r"^Vault",r"^Rn\d",r"^Gm\d+$",r"_rRNA",r"^SNOR",
         r"^Scarna",r"^Rmrp",r"^Rpph",r"^[a-z]?-?R5s"]
import re
def is_noise(g):
    g=str(g); return any(re.search(p,g,re.I) for p in NOISE)

def _trigamma_inverse(x):
    x=np.asarray(x,float); y=0.5+1.0/x
    y=np.where(x>1e7,1.0/np.sqrt(x),y); y=np.where(x<1e-6,1.0/x,y)
    for _ in range(50):
        tri=polygamma(1,y); d=tri*(1-tri/x)/polygamma(2,y); y=y+d
        if np.nanmax(np.abs(d/y))<1e-8: break
    return y

def _ebayes_prior(s2, dfres):
    s2=np.asarray(s2,float); ok=s2>0; z=np.log(s2[ok])
    e=z-digamma(dfres/2)+np.log(dfres/2); ebar=e.mean()
    var_e=np.var(e,ddof=1)-polygamma(1,dfres/2)
    if var_e>0:
        d0=2*_trigamma_inverse(np.array([var_e]))[0]
        s0_2=np.exp(ebar+digamma(d0/2)-np.log(d0/2))
    else:
        d0=np.inf; s0_2=np.exp(ebar)
    return d0,s0_2

def paired_moderated_de(logmat, group_a, group_b, pairs):
    """logmat: genes x samples (log2). pairs: list of (a_sample, b_sample) matched by block.
    Returns DataFrame with log2FC, moderated t, P, FDR. Valid for small n (e.g. n=2)."""
    A=np.column_stack([logmat[a].values for a,_ in pairs])
    B=np.column_stack([logmat[b].values for _,b in pairs])
    D=A-B; n=D.shape[1]
    mean=D.mean(axis=1)
    s2=D.var(axis=1,ddof=1) if n>2 else ((D[:,0]-D[:,1])**2/2)
    dfres=n-1
    d0,s0_2=_ebayes_prior(s2,dfres)
    s2_mod=(d0*s0_2+dfres*s2)/(d0+dfres); df_tot=d0+dfres
    se=np.sqrt(s2_mod/n); t=mean/se
    p=2*stats.t.sf(np.abs(t),df=df_tot)
    order=np.argsort(p); ranks=np.empty_like(order); ranks[order]=np.arange(1,len(p)+1)
    q=p*len(p)/ranks; q=np.minimum.accumulate(q[order][::-1])[::-1]
    fdr=np.empty_like(q); fdr[order]=np.clip(q,0,1)
    return pd.DataFrame({"log2FC":mean,"mod_t":t,"P":p,"FDR":fdr}, index=logmat.index).sort_values("P")

def build_logmatrix(mat, pseudocount=1.0):
    return np.log2(mat+pseudocount)
