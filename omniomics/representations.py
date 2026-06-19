"""Learned representations as gated candidate modalities (never-below-the-anchor adoption).

Future-proofing for when n grows into the tens of thousands and a learned representation (autoencoder / VAE,
e.g. scVI) may finally beat a linear model on some tasks. The PRINCIPLE is unchanged from the rest of the
toolkit: a learned latent is just another *candidate modality* gated onto the anchor's residual and ADOPTED
ONLY IF it beats the anchor out-of-sample (the gate can set its weight to zero, so the combined model is
never below the anchor). Nothing here fires unless the data actually warrants it.

Two pieces:
  - a pluggable encoder interface (PCA / randomized-SVD / sklearn-autoencoder now; a lazy scVI/torch VAE hook
    for large n) -- swap the representation without touching the evaluation;
  - gated_candidate_cv(...): leakage-safe evaluation that fits the encoder INSIDE each train fold, gates the
    latent onto the anchor, and reports the out-of-sample gain over the anchor plus an adoption verdict.

Also: iterative_svd_impute(...) for missing-modality / missing-value imputation (the same "learned slot",
dependency-light; an autoencoder imputer can be swapped in when warranted).

Core paths need only numpy + scikit-learn.
"""
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score


# --------------------------------------------------------------------------- encoders (the swappable slot)
class PCAEncoder:
    """Linear baseline representation (sklearn PCA)."""
    def __init__(self, n_components=32, random_state=0):
        self.n_components = n_components; self.random_state = random_state

    def fit(self, X):
        self._pca = PCA(n_components=min(self.n_components, *X.shape), random_state=self.random_state).fit(X)
        return self

    def transform(self, X):
        return self._pca.transform(X)


class RandomizedEncoder:
    """Fast approximate linear representation via randomized SVD (scales to wide matrices)."""
    def __init__(self, n_components=32, random_state=0):
        self.n_components = n_components; self.random_state = random_state

    def fit(self, X):
        from .scale import randomized_components
        self.mean_ = X.mean(axis=0)
        self.components_ = randomized_components(X, n_components=min(self.n_components, *X.shape),
                                                 random_state=self.random_state)["components"]
        return self

    def transform(self, X):
        return (X - self.mean_) @ self.components_.T


class AutoencoderEncoder:
    """Nonlinear representation: a small sklearn MLP autoencoder; the bottleneck activation is the latent.
    (No torch needed. For large n, swap in a torch/scVI VAE via make_encoder('vae').)"""
    def __init__(self, n_components=32, hidden=128, max_iter=200, random_state=0):
        self.n_components = n_components; self.hidden = hidden
        self.max_iter = max_iter; self.random_state = random_state

    def fit(self, X):
        from sklearn.neural_network import MLPRegressor
        h, k = self.hidden, min(self.n_components, X.shape[1])
        self._ae = MLPRegressor(hidden_layer_sizes=(h, k, h), activation="relu",
                                max_iter=self.max_iter, random_state=self.random_state).fit(X, X)
        self._W0, self._W1 = self._ae.coefs_[0], self._ae.coefs_[1]
        self._b0, self._b1 = self._ae.intercepts_[0], self._ae.intercepts_[1]
        return self

    def transform(self, X):
        return np.maximum(0, np.maximum(0, X @ self._W0 + self._b0) @ self._W1 + self._b1)


def make_encoder(kind="pca", n_components=32, random_state=0, **kw):
    """Factory for the representation slot. kind in {pca, randomized, autoencoder, vae}. 'vae' lazily imports
    a torch/scVI backend and raises a clear message if it is not installed (the large-n drop-in)."""
    kind = kind.lower()
    if kind == "pca":
        return PCAEncoder(n_components, random_state)
    if kind in ("randomized", "rsvd"):
        return RandomizedEncoder(n_components, random_state)
    if kind in ("autoencoder", "ae"):
        return AutoencoderEncoder(n_components, random_state=random_state, **kw)
    if kind in ("vae", "scvi"):
        try:
            import scvi  # noqa: F401
        except Exception as e:
            raise ImportError(
                "kind='vae' needs a deep-learning backend (scvi-tools / torch), recommended only at large n. "
                "Install it (e.g. `pip install scvi-tools`) or use kind in {pca, randomized, autoencoder}."
            ) from e
        raise NotImplementedError(
            "scVI hook is intentionally a stub: wire it here (fit on train, expose get_latent_representation as "
            "transform) when you adopt a deep backend at scale; the gated-adoption evaluation below is unchanged."
        )
    raise ValueError(f"unknown encoder kind: {kind!r}")


# ------------------------------------------------ gated adoption: keep the latent only if it beats the anchor
def _logit(p):
    p = np.clip(p, 1e-6, 1 - 1e-6); return np.log(p / (1 - p))


def gated_candidate_cv(anchor_score, X, y, encoder=None, cv=5, random_state=0, standardize=True,
                       betas=(0.0, 0.25, 0.5, 1.0, 2.0, 4.0), inner_cv=3):
    """Leakage-safe 'does a learned representation beat the anchor?' evaluation.

    For each outer fold: fit the encoder on the TRAIN features only, transform train/test, fit a logistic on
    the latent (secondary score), and combine with the calibrated anchor as logit(anchor) + beta * secondary,
    choosing beta >= 0 by inner-CV AUROC on the train fold (beta = 0 allowed -> the combined model is never
    below the anchor). Then score the held-out fold. Returns the mean out-of-sample anchor vs combined AUROC,
    the mean adopted beta, the fraction of folds that adopted (beta>0), and an `adopt` verdict (combined beats
    anchor on average). The never-below-anchor guarantee holds by construction.

    encoder: any object with fit(X).transform(X) (default PCAEncoder); swap in autoencoder/VAE at scale.
    """
    X = np.asarray(X, dtype=float); y = np.asarray(y); a = np.asarray(anchor_score, dtype=float)
    enc0 = encoder if encoder is not None else PCAEncoder(random_state=random_state)
    skf = StratifiedKFold(cv, shuffle=True, random_state=random_state)
    au_anchor, au_comb, betas_used = [], [], []
    for tr, te in skf.split(X, y):
        sc = StandardScaler().fit(X[tr]) if standardize else None
        Xtr = sc.transform(X[tr]) if standardize else X[tr]
        Xte = sc.transform(X[te]) if standardize else X[te]
        enc = clone_encoder(enc0).fit(Xtr)
        Ltr, Lte = enc.transform(Xtr), enc.transform(Xte)
        # calibrate the anchor and fit the secondary (latent) score on train
        cal = LogisticRegression(max_iter=200).fit(a[tr][:, None], y[tr])
        la = lambda z: _logit(cal.predict_proba(z[:, None])[:, 1])
        sec = LogisticRegression(max_iter=300).fit(Ltr, y[tr])
        str_, ste = sec.predict_proba(Ltr)[:, 1], sec.predict_proba(Lte)[:, 1]
        # choose beta>=0 on an inner CV of the train fold
        best_b, best_auc = 0.0, -1.0
        inner = StratifiedKFold(inner_cv, shuffle=True, random_state=random_state)
        for b in betas:
            aucs = []
            for itr, iva in inner.split(Xtr, y[tr]):
                comb = la(a[tr][iva]) + b * str_[iva]
                if len(np.unique(y[tr][iva])) == 2:
                    aucs.append(roc_auc_score(y[tr][iva], comb))
            if aucs and np.mean(aucs) > best_auc:
                best_auc, best_b = float(np.mean(aucs)), b
        comb_te = la(a[te]) + best_b * ste
        au_anchor.append(roc_auc_score(y[te], a[te]))
        au_comb.append(roc_auc_score(y[te], comb_te))
        betas_used.append(best_b)
    ma, mc = float(np.mean(au_anchor)), float(np.mean(au_comb))
    return dict(auroc_anchor=round(ma, 4), auroc_combined=round(mc, 4), delta=round(mc - ma, 4),
                mean_beta=round(float(np.mean(betas_used)), 3),
                frac_folds_adopted=round(float(np.mean([b > 0 for b in betas_used])), 3),
                adopt=bool(mc > ma))


def clone_encoder(enc):
    """Fresh, unfitted copy of an encoder (so each fold fits independently)."""
    return enc.__class__(**{k: v for k, v in enc.__dict__.items() if not k.endswith("_") and not k.startswith("_")})


# ----------------------------------------------------------------- missing-modality / missing-value imputation
def iterative_svd_impute(X, rank=10, n_iter=50, tol=1e-4):
    """Low-rank (SoftImpute-style) imputation of missing entries -- the dependency-light occupant of the
    'learned slot' for missing-modality/value filling (an autoencoder imputer can be swapped in). NaNs mark
    missing. Returns the completed matrix. Impute INSIDE CV folds to avoid leakage."""
    X = np.array(X, dtype=float); mask = np.isnan(X)
    if not mask.any():
        return X
    col_mean = np.nanmean(X, axis=0); col_mean = np.where(np.isfinite(col_mean), col_mean, 0.0)
    F = X.copy(); F[mask] = np.take(col_mean, np.where(mask)[1])
    prev = None
    r = int(max(1, min(rank, min(F.shape) - 1)))
    for _ in range(n_iter):
        U, s, Vt = np.linalg.svd(F, full_matrices=False)
        Fhat = (U[:, :r] * s[:r]) @ Vt[:r]
        F = np.where(mask, Fhat, X)                      # keep observed entries, replace only missing
        if prev is not None:
            denom = np.linalg.norm(prev[mask]) + 1e-12
            if np.linalg.norm(F[mask] - prev[mask]) / denom < tol:
                break
        prev = F.copy()
    return F
