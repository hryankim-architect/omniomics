"""Unit tests for the out-of-core scale utilities (omniomics.scale). Run in CI; no large data needed."""
import numpy as np
import pytest
from omniomics import scale


def _data(seed=0, n=80, p=600, k_sig=8):
    rng = np.random.default_rng(seed)
    sig = rng.normal(size=n)
    y = (sig + 0.5 * rng.normal(size=n) > 0).astype(int)
    X = rng.normal(size=(n, p))
    X[:, :k_sig] += 1.5 * sig[:, None]          # first k_sig features carry the signal
    return X, y


def test_chunked_columns_reconstructs():
    X, _ = _data()
    blocks = list(scale.chunked_columns(X, chunk_size=128))
    cols = np.concatenate([idx for idx, _ in blocks])
    recon = np.concatenate([B for _, B in blocks], axis=1)
    assert np.array_equal(cols, np.arange(X.shape[1]))
    assert np.allclose(recon, X)


def test_chunked_columns_npy_roundtrip(tmp_path):
    X, _ = _data()
    p = tmp_path / "X.npy"; np.save(p, X)
    recon = np.concatenate([B for _, B in scale.chunked_columns(str(p), chunk_size=100)], axis=1)
    assert np.allclose(recon, X)


def test_streaming_sis_recovers_signal_and_matches_inmemory():
    X, y = _data()
    out = scale.streaming_sis(X, y, top_k=20, chunk_size=64)
    # the 8 signal features must be in the top-20
    assert len(set(range(8)) & set(out["indices"].tolist())) == 8
    # streaming score == in-memory score (residualized-correlation), same ranking
    yc = y - y.mean()
    Xc = X - X.mean(0)
    full = np.abs(Xc.T @ yc) / (np.sqrt((Xc * Xc).sum(0)) * np.sqrt(yc @ yc) + 1e-12)
    top_mem = np.sort(np.argsort(-full)[:20])
    assert np.array_equal(out["indices"], top_mem)


def test_streaming_gram_pca_matches_exact_pca():
    X, _ = _data(seed=1, n=60, p=300)
    k = 5
    g = scale.streaming_gram_pca(X, n_components=k, chunk_size=64)
    # compare against exact SVD sample scores (sign-agnostic): |corr| per component ~ 1
    Xc = X - X.mean(0)
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    exact = (U * S)[:, :k]
    for j in range(k):
        c = abs(np.corrcoef(g["scores"][:, j], exact[:, j])[0, 1])
        assert c > 0.999
    assert np.allclose(np.sort(g["singular_values"])[::-1], np.sort(S[:k])[::-1], rtol=1e-6, atol=1e-6)


def test_randomized_components_approximates_svd():
    # low-rank-plus-noise (realistic omics-like spectrum with a clear gap), where randomized SVD is accurate
    rng = np.random.default_rng(2); n, p, r0 = 80, 400, 6
    F = rng.normal(size=(n, r0)) * np.array([6.0, 5.0, 4.0, 3.0, 2.5, 2.0])
    L = rng.normal(size=(r0, p))
    X = F @ L + 0.3 * rng.normal(size=(n, p))
    k = r0
    out = scale.randomized_components(X, n_components=k, random_state=0)
    Xc = X - X.mean(0)
    _, S, _ = np.linalg.svd(Xc, full_matrices=False)
    assert np.allclose(out["singular_values"], S[:k], rtol=0.02)      # leading k well-separated -> accurate
    assert out["scores"].shape == (n, k) and out["components"].shape == (k, p)
