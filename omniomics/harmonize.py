"""Cross-study harmonization — the Phase 1 core.
Combine heterogeneous per-study expression matrices (different units, gene IDs,
pipelines) onto a common gene space, then correct study batch effects so biology,
not provenance, drives downstream structure."""
import numpy as np, pandas as pd

def common_gene_panel(study_mats: dict[str, pd.DataFrame], min_studies=None):
    """study_mats: {study_id: genes(symbol,UPPER) x samples}. Intersect genes present in all (or >=min_studies)."""
    from collections import Counter
    c = Counter()
    for m in study_mats.values():
        c.update(set(m.index))
    need = len(study_mats) if min_studies is None else min_studies
    genes = sorted([g for g, k in c.items() if k >= need])
    return genes

def assemble(study_mats: dict[str, pd.DataFrame], genes=None, log=True):
    """Return combined matrix (genes x all-samples) + sample->study map."""
    if genes is None:
        genes = common_gene_panel(study_mats)
    cols, study_of = [], {}
    blocks = []
    for sid, m in study_mats.items():
        sub = m.reindex(genes)
        sub = sub.rename(columns={c: f"{sid}|{c}" for c in sub.columns})
        for c in sub.columns: study_of[c] = sid
        blocks.append(sub)
    X = pd.concat(blocks, axis=1)
    if log:
        X = np.log2(X.clip(lower=0) + 1)
    return X, pd.Series(study_of)

def quantile_normalize(X: pd.DataFrame):
    """Make sample distributions comparable across studies/units."""
    ranks = X.rank(method="average")
    mean_sorted = np.sort(X.values, axis=0).mean(axis=1)
    qn = ranks.apply(lambda col: np.interp(col, np.arange(1, len(col)+1), mean_sorted))
    return pd.DataFrame(qn.values, index=X.index, columns=X.columns)

def combat_lite(X: pd.DataFrame, study_of: pd.Series):
    """Location/scale batch correction (ComBat-style, no EB shrinkage): per-study
    center each gene to the grand mean and rescale to pooled variance."""
    Xc = X.copy()
    grand_mean = X.mean(axis=1)
    pooled_sd = X.std(axis=1).replace(0, 1)
    for sid in study_of.unique():
        cols = study_of[study_of == sid].index
        b = X[cols]
        bm = b.mean(axis=1); bsd = b.std(axis=1).replace(0, 1)
        Xc[cols] = ((b.sub(bm, axis=0)).div(bsd, axis=0)).mul(pooled_sd, axis=0).add(grand_mean, axis=0)
    return Xc

def combat_eb(X: pd.DataFrame, batch, mod=None, max_iter=50):
    """Parametric empirical-Bayes ComBat (Johnson, Li & Rabinovic 2007).
    Unlike combat_lite (location/scale per batch), EB shrinks per-batch effects toward a
    common prior — the right tool for N>2 batches and small/uneven batches. `mod` (samples x p,
    no intercept) is a biological covariate to PRESERVE (e.g., subtype), which is what makes
    correction non-trivial when batch is confounded with biology."""
    genes, samples = X.index, X.columns
    Y = X.values.astype(float)                         # g x n
    batch = np.asarray(batch); batches = list(pd.unique(batch))
    idxs = [np.where(batch == b)[0] for b in batches]; nb = len(batches); n = Y.shape[1]
    Xb = np.zeros((n, nb))
    for j, ix in enumerate(idxs): Xb[ix, j] = 1
    design = Xb if mod is None else np.hstack([Xb, np.asarray(mod, float)])
    B = np.linalg.lstsq(design, Y.T, rcond=None)[0]    # p x g
    grand = (np.array([len(ix) for ix in idxs]) / n) @ B[:nb]      # g
    resid = Y.T - design @ B
    var_pooled = (resid ** 2).mean(axis=0)             # g
    stand = np.tile(grand[:, None], (1, n))
    if mod is not None: stand = stand + (np.asarray(mod, float) @ B[nb:]).T
    sd = np.sqrt(var_pooled + 1e-12)[:, None]
    Z = (Y - stand) / sd
    g_hat = np.vstack([Z[:, ix].mean(axis=1) for ix in idxs])
    d_hat = np.vstack([Z[:, ix].var(axis=1, ddof=1) + 1e-8 for ix in idxs])
    g_bar = g_hat.mean(axis=1); t2 = g_hat.var(axis=1, ddof=1)
    def aprior(d): m = d.mean(); s2 = d.var(); return (2 * s2 + m * m) / s2
    def bprior(d): m = d.mean(); s2 = d.var(); return (m * s2 + m ** 3) / s2
    Zadj = Z.copy()
    for j, ix in enumerate(idxs):
        a, b = aprior(d_hat[j]), bprior(d_hat[j]); ni = len(ix)
        g, d = g_hat[j].copy(), d_hat[j].copy()
        for _ in range(max_iter):
            gn = (t2[j] * ni * g_hat[j] + d * g_bar[j]) / (t2[j] * ni + d)
            s2 = ((Z[:, ix] - gn[:, None]) ** 2).sum(axis=1)
            dn = (0.5 * s2 + b) / (ni / 2.0 + a - 1)
            if np.max(np.abs(gn - g) / (np.abs(g) + 1e-8)) < 1e-4 and \
               np.max(np.abs(dn - d) / (np.abs(d) + 1e-8)) < 1e-4:
                g, d = gn, dn; break
            g, d = gn, dn
        Zadj[:, ix] = (Z[:, ix] - g[:, None]) / np.sqrt(d)[:, None]
    Yadj = Zadj * sd + stand
    return pd.DataFrame(Yadj, index=genes, columns=samples)

def batch_pred_accuracy(X: pd.DataFrame, batch):
    """How recoverable is the batch label from the data (lower = better mixing)."""
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import cross_val_score
    return cross_val_score(RandomForestClassifier(n_estimators=120, random_state=0),
                           X.T.values, np.asarray(batch), cv=4).mean()

def pca(X: pd.DataFrame, n=10):
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    Z = StandardScaler().fit_transform(X.T.values)         # samples x genes
    p = PCA(n_components=min(n, Z.shape[0]-1)).fit(Z)
    return p.transform(Z), p.explained_variance_ratio_

def batch_variance_explained(X: pd.DataFrame, study_of: pd.Series):
    """Fraction of PC1-2 variance attributable to study label (a batch-effect score).
    Lower after correction = better mixing."""
    pcs, evr = pca(X, n=2)
    df = pd.DataFrame(pcs, columns=["PC1","PC2"]); df["study"] = study_of.values
    out = {}
    for pc in ["PC1","PC2"]:
        grand = df[pc].mean()
        ss_tot = ((df[pc]-grand)**2).sum()
        ss_between = sum(len(g)*(g[pc].mean()-grand)**2 for _,g in df.groupby("study"))
        out[pc] = ss_between/ss_tot if ss_tot>0 else 0.0
    return out, evr
