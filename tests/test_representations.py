"""Unit tests for learned-representation gated adoption (omniomics.representations). CI-fast, no torch."""
import numpy as np
import pytest
from omniomics import representations as rep


def test_encoders_shapes():
    rng = np.random.default_rng(0); X = rng.normal(size=(60, 50))
    for enc in (rep.PCAEncoder(8), rep.RandomizedEncoder(8), rep.AutoencoderEncoder(8, hidden=32, max_iter=80)):
        L = enc.fit(X).transform(X)
        assert L.shape == (60, 8)


def test_make_encoder_and_vae_hook():
    assert isinstance(rep.make_encoder("pca", 4), rep.PCAEncoder)
    assert isinstance(rep.make_encoder("randomized", 4), rep.RandomizedEncoder)
    assert isinstance(rep.make_encoder("autoencoder", 4), rep.AutoencoderEncoder)
    with pytest.raises((ImportError, NotImplementedError)):   # no torch/scvi installed -> clear failure
        rep.make_encoder("vae", 4)


def test_never_below_anchor_when_gate_zero():
    # with betas=(0.0,) the combined score is a monotone transform of the anchor -> identical AUROC (the
    # structural never-below-anchor guarantee).
    rng = np.random.default_rng(1); n = 200
    sig = rng.normal(size=n); y = (sig + 0.5 * rng.normal(size=n) > 0).astype(int)
    anchor = sig + 0.3 * rng.normal(size=n)
    X = rng.normal(size=(n, 30))
    out = rep.gated_candidate_cv(anchor, X, y, encoder=rep.PCAEncoder(8), cv=4, betas=(0.0,), inner_cv=2)
    assert out["auroc_combined"] == out["auroc_anchor"]
    assert out["adopt"] is False and out["frac_folds_adopted"] == 0.0


def test_adopts_only_when_representation_helps():
    rng = np.random.default_rng(2); n = 260
    sig = rng.normal(size=n); y = (sig > 0).astype(int)
    anchor = 0.4 * sig + rng.normal(size=n)                 # weak, noisy anchor
    X = rng.normal(size=(n, 40)); X[:, :6] += 2.0 * sig[:, None]   # X carries the clean signal the anchor lacks
    helped = rep.gated_candidate_cv(anchor, X, y, encoder=rep.PCAEncoder(8), cv=4, inner_cv=2, random_state=0)
    assert helped["adopt"] and helped["delta"] > 0 and helped["frac_folds_adopted"] > 0

    # null case: X is pure noise -> the gate should not adopt a useless representation (delta ~ 0)
    Xn = rng.normal(size=(n, 40))
    null = rep.gated_candidate_cv(anchor, Xn, y, encoder=rep.PCAEncoder(8), cv=4, inner_cv=2, random_state=0)
    assert null["delta"] <= helped["delta"] and abs(null["delta"]) < 0.05


def test_iterative_svd_impute_recovers_low_rank():
    rng = np.random.default_rng(3); n, p, r = 60, 40, 3
    truth = rng.normal(size=(n, r)) @ rng.normal(size=(r, p))
    Xo = truth.copy(); mask = rng.random((n, p)) < 0.2; Xo[mask] = np.nan
    filled = rep.iterative_svd_impute(Xo, rank=r, n_iter=80)
    assert not np.isnan(filled).any()
    rmse = np.sqrt(np.mean((filled[mask] - truth[mask]) ** 2))
    base = np.sqrt(np.mean((np.nanmean(truth) - truth[mask]) ** 2))   # mean-fill baseline
    assert rmse < 0.5 * base                                          # low-rank impute clearly beats mean-fill
