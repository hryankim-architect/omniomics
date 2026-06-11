"""Joint multi-omics integration: MOFA-style linear factor model + pathway scoring.
For continuous omics (log-expression, methylation beta) on matched samples — the regime
where MOFA/iCluster apply (scVI/PeakVI need raw counts, so are not used here)."""
import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def _scale_impute(df):
    X = df.T.copy()                                  # samples x genes
    X = X.fillna(X.mean())
    Z = StandardScaler().fit_transform(X.values)
    return Z                                          # samples x genes

def mofa_lite(views: dict, k=10):
    """views: {name: genes x samples (matched samples, same column order)}.
    Returns factor scores (samples x k), and per-view variance-explained table (view x factor)."""
    names = list(views)
    samples = views[names[0]].columns
    blocks, spans = [], {}
    start = 0
    for nm in names:
        Z = _scale_impute(views[nm])                 # samples x genes_v
        spans[nm] = (start, start + Z.shape[1]); start += Z.shape[1]
        blocks.append(Z)
    C = np.hstack(blocks)                             # samples x all-features
    p = PCA(n_components=k).fit(C)
    T = p.transform(C)                               # samples x k (factor scores)
    P = p.components_.T                               # features x k (loadings)
    lam = T.var(axis=0)                              # variance of each factor
    r2 = {}
    for nm in names:
        s, e = spans[nm]; nfeat = e - s
        load_sq = (P[s:e, :] ** 2).sum(axis=0)       # per-factor loading energy in view
        r2[nm] = lam * load_sq / nfeat               # fraction of view variance per factor
    R2 = pd.DataFrame(r2, index=[f"F{i+1}" for i in range(k)]).T
    Tdf = pd.DataFrame(T, index=samples, columns=[f"F{i+1}" for i in range(k)])
    return Tdf, R2

def dmoi_representation(rna, meth, poles: dict, gmt: dict, meth_sign=-1.0):
    """DMOI-style structured fusion (Swanson/Kim dmoi-brca design, generalized).
    For each pole, build a pole-conditioned RNA score, a methylation score (meth_sign=-1 so
    high methylation -> low-expression direction), and an RNA-vs-methylation DISAGREEMENT scalar.
    rna, meth: genes(UPPER) x samples (matched columns). poles: {pole_name: [HALLMARK_set,...]}.
    Returns samples x features DataFrame with columns grouped rna_/meth_/disagree_."""
    from sklearn.preprocessing import StandardScaler
    S = rna.columns
    zR = pd.DataFrame(StandardScaler().fit_transform(rna.T.fillna(rna.T.mean()).values), index=S, columns=rna.index)
    zM = pd.DataFrame(StandardScaler().fit_transform((meth_sign*meth).T.fillna((meth_sign*meth).T.mean()).values), index=S, columns=meth.index)
    cols = {}
    for p, sets in poles.items():
        genes = set()
        for s in sets: genes |= set(gmt.get(s, []))
        gr = [g for g in genes if g in rna.index]; gm = [g for g in genes if g in meth.index]
        r = zR[gr].mean(axis=1) if gr else pd.Series(0, index=S)
        m = zM[gm].mean(axis=1) if gm else pd.Series(0, index=S)
        cols[f"rna_{p}"] = r
        cols[f"meth_{p}"] = m
        cols[f"disagree_{p}"] = (r - m)            # where RNA and methylation diverge
    return pd.DataFrame(cols, index=S)

def load_gmt(path):
    sets = {}
    for L in open(path):
        parts = L.rstrip("\n").split("\t")
        if len(parts) >= 3:
            sets[parts[0]] = [g.upper() for g in parts[2:] if g]
    return sets

def pathway_scores(view: pd.DataFrame, gmt: dict, min_genes=5):
    """view: genes(UPPER) x samples. Returns pathways x samples = mean z-score of member genes."""
    Z = pd.DataFrame(StandardScaler().fit_transform(view.T.fillna(view.T.mean()).values),
                     index=view.columns, columns=view.index)   # samples x genes
    rows = {}
    for pw, genes in gmt.items():
        g = [x for x in genes if x in view.index]
        if len(g) >= min_genes:
            rows[pw] = Z[g].mean(axis=1)
    return pd.DataFrame(rows).T                         # pathways x samples
