"""Out-of-core loading + randomized linear algebra for genome-wide (p >> n) scale.

Future-proofing utilities so the anchored-discovery toolkit can run on feature matrices that do not fit in
memory (e.g. EPIC methylation, ~850k probes) without changing the science:

  - chunked_columns(source, ...)   : stream a feature matrix in column blocks from an ndarray / DataFrame /
                                     .npy (mmap) / .parquet / HDF5 -- one block in memory at a time.
  - streaming_sis(source, y, ...)  : one-pass, out-of-core Sure-Independence-Screening (rank features by
                                     |corr with the anchor-residualized response|) -- the disk-backed twin
                                     of the in-memory screen in multiomics.anchored_residual_discovery.
  - streaming_gram_pca(source, ...): out-of-core PCA for p >> n via the n x n Gram matrix X Xᵀ accumulated
                                     over column blocks (never holds all p features), then eigendecomposed.
  - randomized_components(X, k)    : thin wrapper over sklearn randomized_svd -- fast approximate top-k.

Core paths need only numpy + scikit-learn. Parquet needs `pyarrow`, HDF5 needs `h5py` (imported lazily with
a clear message). Designed to be drop-in for screen_top: screen out-of-core, then load only the survivors.
"""
import numpy as np
import pandas as pd


def _ncols(source):
    if isinstance(source, np.ndarray):
        return source.shape[1]
    if isinstance(source, pd.DataFrame):
        return source.shape[1]
    if isinstance(source, tuple):                      # (hdf5_path, dataset)
        import h5py  # noqa
        with __import__("h5py").File(source[0], "r") as f:
            return f[source[1]].shape[1]
    s = str(source)
    if s.endswith(".npy"):
        return np.load(s, mmap_mode="r").shape[1]
    if s.endswith(".parquet"):
        import pyarrow.parquet as pq
        return pq.ParquetFile(s).schema_arrow.num_fields
    raise ValueError(f"unsupported source for chunked loading: {source!r}")


def feature_names(source):
    """Best-effort feature (column) names; falls back to integer indices."""
    if isinstance(source, pd.DataFrame):
        return list(source.columns)
    s = str(source)
    if s.endswith(".parquet"):
        import pyarrow.parquet as pq
        return list(pq.ParquetFile(s).schema_arrow.names)
    return list(range(_ncols(source)))


def chunked_columns(source, chunk_size=2048):
    """Yield (col_indices ndarray, block ndarray of shape (n_samples, n_cols_in_block)).

    source: np.ndarray (samples x features), pd.DataFrame, path to .npy (memory-mapped), path to .parquet
    (needs pyarrow), or (hdf5_path, dataset_name) (needs h5py). Only one block is materialised at a time.
    """
    p = _ncols(source)
    s = str(source)
    if isinstance(source, np.ndarray):
        get = lambda lo, hi: np.asarray(source[:, lo:hi], dtype=float)
    elif isinstance(source, pd.DataFrame):
        get = lambda lo, hi: source.iloc[:, lo:hi].to_numpy(dtype=float)
    elif isinstance(source, tuple):
        import h5py
        f = h5py.File(source[0], "r"); dset = f[source[1]]
        get = lambda lo, hi: np.asarray(dset[:, lo:hi], dtype=float)
    elif s.endswith(".npy"):
        mm = np.load(s, mmap_mode="r")
        get = lambda lo, hi: np.asarray(mm[:, lo:hi], dtype=float)
    elif s.endswith(".parquet"):
        import pyarrow.parquet as pq
        pf = pq.ParquetFile(s); names = pf.schema_arrow.names
        get = lambda lo, hi: pf.read(columns=names[lo:hi]).to_pandas().to_numpy(dtype=float)
    else:
        raise ValueError(f"unsupported source: {source!r}")
    for lo in range(0, p, chunk_size):
        hi = min(lo + chunk_size, p)
        yield np.arange(lo, hi), get(lo, hi)


def _residualize(y, anchor):
    """Return the anchor-residualized, mean-centered response (or just centered y if anchor is None)."""
    y = np.asarray(y, dtype=float); yc = y - y.mean()
    if anchor is None:
        return yc
    a = np.asarray(anchor, dtype=float); ac = a - a.mean()
    den = float(ac @ ac) + 1e-12
    r = yc - (float(ac @ yc) / den) * ac
    return r - r.mean()


def streaming_sis(source, y, anchor=None, top_k=5000, chunk_size=2048):
    """One-pass out-of-core Sure-Independence-Screening.

    Rank every feature by |Pearson correlation with the anchor-residualized response| while streaming the
    matrix in column blocks, and return the top_k. Lets SIS scale to feature matrices that do not fit in
    memory; downstream you load only the survivors. Returns dict(indices, scores, names).
    """
    rY = _residualize(y, anchor); rYn = np.sqrt(float(rY @ rY)) + 1e-12
    p = _ncols(source); scores = np.empty(p, dtype=float)
    for idx, B in chunked_columns(source, chunk_size):
        Bc = B - B.mean(axis=0)
        num = np.abs(Bc.T @ rY)
        den = np.sqrt((Bc * Bc).sum(axis=0)) * rYn + 1e-12
        scores[idx] = num / den
    k = min(int(top_k), p)
    keep = np.sort(np.argsort(-scores)[:k])
    names = feature_names(source)
    return dict(indices=keep, scores=scores[keep], names=[names[i] for i in keep])


def streaming_gram_pca(source, n_components=50, chunk_size=2048, center=True):
    """Out-of-core PCA for p >> n: accumulate the n x n Gram matrix X Xᵀ over column blocks (never holding
    all p features in memory), then eigendecompose. Returns dict(scores (n x k), explained_variance (k,),
    singular_values (k,)). Feature-centering (each feature mean 0 across samples) uses one extra streaming
    pass over the column means.
    """
    # pass 1: feature (column) means
    if center:
        p = _ncols(source); mu = np.empty(p, dtype=float)
        for idx, B in chunked_columns(source, chunk_size):
            mu[idx] = B.mean(axis=0)
    # pass 2: Gram of centered blocks
    G = None
    for idx, B in chunked_columns(source, chunk_size):
        Bc = B - mu[idx] if center else B
        G = (Bc @ Bc.T) if G is None else G + (Bc @ Bc.T)
    n = G.shape[0]
    w, V = np.linalg.eigh(G)                            # ascending
    order = np.argsort(w)[::-1][:int(n_components)]
    w = np.clip(w[order], 0, None); V = V[:, order]
    sv = np.sqrt(w)
    scores = V * sv                                     # U * S  (n x k sample embedding)
    return dict(scores=scores, explained_variance=w / max(n - 1, 1), singular_values=sv)


def randomized_components(X, n_components=50, random_state=0, n_oversamples=10, n_iter=7):
    """Thin wrapper over sklearn randomized_svd: fast approximate top-k SVD of an in-memory (n x p) matrix
    (feature-centered). Returns dict(scores (n x k) = U*S, components (k x p) = Vᵀ, singular_values (k,)).
    For matrices that do not fit in memory, use streaming_gram_pca instead.
    """
    from sklearn.utils.extmath import randomized_svd
    X = np.asarray(X, dtype=float); Xc = X - X.mean(axis=0)
    U, S, Vt = randomized_svd(Xc, n_components=int(n_components), n_oversamples=n_oversamples,
                              n_iter=n_iter, random_state=random_state)
    return dict(scores=U * S, components=Vt, singular_values=S)
