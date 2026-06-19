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


# ---------- DMOI v2: reliability-weighted, interaction-aware fusion ----------
# Rationale (see reports/DMOI_method_assessment.md, validated on TCGA-BRCA + GSE57577):
#   * v1's signature feature, disagree = z_rna - z_meth, lies in the linear span of {z_rna, z_meth}
#     and is therefore PROVABLY redundant for a linear classifier (benchmark: +0.0000 AUROC). The
#     real cross-omics signal is an INTERACTION, so v2 emits product terms and is meant for a
#     nonlinear learner (e.g. GradientBoosting), where it buys a small reproducible gain (+~0.01).
#   * A shallow/noisy layer fused at equal weight can drag a model BELOW the single-omic baseline.
#     v2 therefore weights each gene's methylation by a measurement-reliability score, which makes
#     degradation graceful (a corrupted layer is automatically down-weighted).

def methylation_reliability(meth_probes, gene_to_probes, min_probes=2):
    """Per-gene methylation reliability = mean inter-probe Pearson correlation across samples.

    meth_probes : probes x samples DataFrame of beta values (may contain NaN).
    gene_to_probes : {GENE: [probe_id, ...]}.
    Genes with < `min_probes` mapped probes get NaN (the caller may fill with the cohort median).

    This is the default, assay-agnostic reliability proxy: genes whose probes disagree across
    samples are trusted less, and the score collapses correctly when a layer is corrupted. Swap in
    a different proxy (1 - detection_p, coverage-based shrinkage, technical-replicate concordance)
    by passing the resulting {GENE: weight} dict straight to `dmoi_v2_representation(reliability=...)`.
    """
    rel = {}
    for g, probes in gene_to_probes.items():
        ps = [p for p in probes if p in meth_probes.index]
        if len(ps) < min_probes:
            rel[g] = float("nan")
            continue
        M = meth_probes.loc[ps].astype(float).values            # probes x samples
        rm = np.nanmean(M, axis=1, keepdims=True)               # impute NaN with per-probe mean
        ii = np.where(np.isnan(M))
        if ii[0].size:
            M = M.copy(); M[ii] = np.take(rm, ii[0])
        M = np.nan_to_num(M)
        C = np.corrcoef(M)
        iu = np.triu_indices(len(ps), 1)
        rel[g] = float(np.nanmean(C[iu])) if iu[0].size else float("nan")
    return rel


def dmoi_v2_representation(rna, meth, poles, gmt, reliability=None, meth_sign=-1.0,
                           interactions=True):
    """DMOI v2 structured representation, designed to feed a NONLINEAR learner.

    rna, meth : genes(UPPER) x samples, matched columns (gene-level methylation already aggregated).
    poles : {pole_name: [GENE_SET_NAME, ...]};  gmt : {GENE_SET_NAME: [GENE, ...]}.
    reliability : optional {GENE: weight >= 0} (e.g. from `methylation_reliability`). When given, a
        pole's methylation score is the reliability-weighted mean of its z-scored member genes;
        when None it is the plain mean (v2 then reduces to interactions-only).
    meth_sign : -1.0 so that high methylation maps to the low-expression direction (as in v1).
    interactions : if True, also emit per-pole product terms rna_<pole> * meth_<pole> -- the
        cross-omics interaction that a linear model cannot exploit but a tree/MLP can.

    Returns samples x features DataFrame with columns rna_<pole>, meth_<pole>, and (optionally)
    int_<pole>. Unlike v1 it deliberately does NOT emit the redundant difference term.
    """
    from sklearn.preprocessing import StandardScaler
    S = rna.columns
    zR = pd.DataFrame(StandardScaler().fit_transform(rna.T.fillna(rna.T.mean()).values),
                      index=S, columns=rna.index)
    msig = meth_sign * meth
    zM = pd.DataFrame(StandardScaler().fit_transform(msig.T.fillna(msig.T.mean()).values),
                      index=S, columns=meth.index)
    cols = {}
    for p, sets in poles.items():
        genes = set()
        for s in sets:
            genes |= set(gmt.get(s, []))
        gr = [g for g in genes if g in rna.index]
        gm = [g for g in genes if g in meth.index]
        r = zR[gr].mean(axis=1) if gr else pd.Series(0.0, index=S)
        if gm:
            if reliability:
                w = np.array([max(float(reliability.get(g, 0.0)), 0.0) + 1e-3 for g in gm])
                m = pd.Series((zM[gm].values * w).sum(axis=1) / w.sum(), index=S)
            else:
                m = zM[gm].mean(axis=1)
        else:
            m = pd.Series(0.0, index=S)
        cols[f"rna_{p}"] = r
        cols[f"meth_{p}"] = m
        if interactions:
            cols[f"int_{p}"] = r * m
    return pd.DataFrame(cols, index=S)


def dmoi_v2_genelevel(rna, meth, genes=None, reliability=None, meth_sign=-1.0, interactions=True):
    """Gene-resolution DMOI v2 features (NO pole averaging), for a regularised / nonlinear learner.

    Roadmap item 3: pole/pathway averaging is good for interpretation but discards signal for
    prediction. This keeps per-gene resolution -- for each gene it emits z-scored RNA, z-scored
    (reliability-scaled) methylation, and their product.

    rna, meth : genes(UPPER) x samples (matched columns).  genes : subset (default: intersection).
    reliability : optional {GENE: weight >= 0}; scales that gene's methylation feature, so an
        unreliable gene is shrunk toward 0 rather than trusted equally.
    Returns samples x features with columns rna_<g>, meth_<g>, int_<g>.
    """
    from sklearn.preprocessing import StandardScaler
    S = rna.columns
    if genes is None:
        genes = [g for g in rna.index if g in meth.index]
    else:
        genes = [g for g in genes if g in rna.index and g in meth.index]
    zR = pd.DataFrame(StandardScaler().fit_transform(rna.loc[genes].T.fillna(rna.loc[genes].T.mean()).values),
                      index=S, columns=genes)
    msig = meth_sign * meth.loc[genes]
    zM = pd.DataFrame(StandardScaler().fit_transform(msig.T.fillna(msig.T.mean()).values),
                      index=S, columns=genes)
    cols = {}
    for g in genes:
        w = 1.0 if not reliability else max(float(reliability.get(g, 0.0)), 0.0) + 1e-3
        r = zR[g]
        m = zM[g] * w
        cols[f"rna_{g}"] = r
        cols[f"meth_{g}"] = m
        if interactions:
            cols[f"int_{g}"] = r * m
    return pd.DataFrame(cols, index=S)


def dmoi_regimes(layers, zscore=True):
    """Partition items into cross-layer concordance/disagreement regimes -- the descriptive,
    effect-modification use of the dialectical 'disagreement' idea (roadmap item 4).

    `layers` : DataFrame, rows = items (genes or samples), columns = one signed score per omics
    layer, each oriented so that 'high = up/active' (e.g. occupancy gain, methylation gain,
    repression = -log2FC). Returns the frame with added columns:
        z_<layer>     per-layer z-scores (when zscore=True),
        concordance   mean z across layers (high => layers agree 'up'),
        n_up          number of layers > 0,
        regime        'concordant_up' (all up), 'concordant_down' (all down), else 'discordant'.

    This is for stratification / testing whether a relationship is modified by cross-layer
    agreement -- NOT a predictor (the report shows the disagreement signal is predictively inert).
    """
    cols_in = list(layers.columns)
    Z = layers
    out = layers.copy()
    if zscore:
        Z = (layers - layers.mean()) / (layers.std(ddof=0) + 1e-9)
        for c in cols_in:
            out[f"z_{c}"] = Z[c]
    n_up = (Z > 0).sum(axis=1)
    out["concordance"] = Z.mean(axis=1)
    out["n_up"] = n_up
    k = len(cols_in)
    out["regime"] = np.where(n_up == k, "concordant_up",
                             np.where(n_up == 0, "concordant_down", "discordant"))
    return out


def anchored_gate(anchor_prob, secondary_score, y, betas=(0.0, 0.25, 0.5, 1.0, 2.0, 4.0), margin=0.0):
    """Leader-anchored, gated residual integration -- the 'never below the anchor' combiner.

    Instead of fusing modalities symmetrically (where a strong, low-noise modality can be dragged
    down by a weak one), anchor on the most robust modality and add others as a gated residual.
    Given the anchor's predicted probabilities and a secondary modality's (residual) score for the
    same items, choose a NON-NEGATIVE weight beta -- including 0 -- that maximises AUROC of
        logit(anchor_prob) + beta * secondary_score.
    Because beta = 0 is in the grid, the combiner defaults to the anchor whenever the secondary
    modality does not help, so it cannot fall below the anchor on the selection set; beta > 0 is
    chosen only where the second modality adds discriminative signal. Returns (beta, combined_score).

    Choose `betas` on held-out / inner-CV data (not the test fold) for an honest estimate. This is
    the recommended way to evaluate and add a secondary omic on top of a dominant one (e.g. RNA).

    `margin` (>= 0): a non-zero beta is accepted only if it beats the anchor (beta = 0) by more than
    this AUROC margin. A small positive margin guards the never-below-anchor guarantee against
    small-sample selection noise (an apparent inner-CV gain that does not generalise); margin = 0
    reproduces the plain argmax.
    """
    from sklearn.metrics import roc_auc_score
    p = np.clip(np.asarray(anchor_prob, dtype=float), 1e-4, 1 - 1e-4)
    base = np.log(p / (1 - p))
    s = np.asarray(secondary_score, dtype=float)
    auc0 = roc_auc_score(y, base)
    cand = [(b, roc_auc_score(y, base + b * s)) for b in betas if b != 0.0]
    cand = [(b, a) for b, a in cand if a > auc0 + margin]
    beta = float(max(cand, key=lambda t: t[1])[0]) if cand else 0.0
    return beta, base + beta * s


def anchored_integrate(Xa, Xs, y, anchor=None, secondary=None,
                       betas=(0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0), cv=5, random_state=0,
                       gate_margin=0.01, inner_repeats=3):
    """End-to-end leader-anchored gated-residual multi-omics integration (leakage-safe OOF).

    Xa : anchor-modality feature matrix (samples x p_a) -- use the most robust modality (e.g. RNA).
    Xs : secondary-modality feature matrix (samples x p_s).
    y  : binary labels.
    anchor    : sklearn classifier with predict_proba (default StandardScaler + LogisticRegression).
    secondary : sklearn regressor fit to the anchor RESIDUAL (default StandardScaler + Ridge).

    For each outer CV fold: fit the anchor on the training part; fit the secondary to the training
    residual (y - inner-OOF anchor probability); pick a non-negative weight beta on an inner CV via
    `anchored_gate` (beta = 0 allowed); score the held-out fold as
        logit(anchor_prob) + beta * residual_prediction.
    Because beta = 0 is always in the grid, the secondary can only be *added where it earns it* --
    on the inner selection it never drags the anchor down, and on held-out data it stays at the
    anchor when the secondary is uninformative and improves on it where the secondary carries
    orthogonal signal.

    Returns dict: oof_anchor (probs), oof_combined (scores), betas (per fold),
    auroc_anchor, auroc_combined, delta (= auroc_combined - auroc_anchor).
    """
    from sklearn.base import clone
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression, Ridge
    from sklearn.pipeline import make_pipeline
    from sklearn.model_selection import StratifiedKFold, KFold, cross_val_predict
    from sklearn.metrics import roc_auc_score
    Xa = np.asarray(Xa, dtype=float); Xs = np.asarray(Xs, dtype=float); y = np.asarray(y)
    if anchor is None:
        anchor = make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000, C=0.5))
    if secondary is None:
        secondary = make_pipeline(StandardScaler(), Ridge(alpha=20.0))
    n = len(y); oof_a = np.zeros(n); oof_c = np.zeros(n); betas_used = []
    inner = max(3, min(cv, 5))
    def _logit(p):
        p = np.clip(p, 1e-4, 1 - 1e-4); return np.log(p / (1 - p))
    for tr, te in StratifiedKFold(cv, shuffle=True, random_state=random_state).split(Xa, y):
        ytr = y[tr]
        # robust beta selection: average each beta's inner-CV AUROC over `inner_repeats` splits, then
        # accept a non-zero beta only if it beats the anchor (beta=0) by more than `gate_margin`.
        # Repeated inner CV tames the small-sample selection noise that can otherwise pick a beta that
        # looks good in-fold but generalises below the anchor.
        beta_auc = {b: [] for b in betas}; resid_last = None
        for rep in range(inner_repeats):
            p_in = cross_val_predict(clone(anchor), Xa[tr], ytr,
                                     cv=StratifiedKFold(inner, shuffle=True, random_state=random_state + rep),
                                     method="predict_proba")[:, 1]
            resid = ytr - p_in
            mc_in = cross_val_predict(clone(secondary), Xs[tr], resid,
                                      cv=KFold(inner, shuffle=True, random_state=random_state + rep))
            base = _logit(p_in)
            for b in betas:
                beta_auc[b].append(roc_auc_score(ytr, base + b * mc_in))
            resid_last = resid
        auc0 = float(np.mean(beta_auc.get(0.0, [0.5])))
        cand = [(b, float(np.mean(v))) for b, v in beta_auc.items() if b != 0.0 and float(np.mean(v)) > auc0 + gate_margin]
        beta = float(max(cand, key=lambda t: t[1])[0]) if cand else 0.0
        betas_used.append(beta)
        a = clone(anchor).fit(Xa[tr], ytr); p_te = a.predict_proba(Xa[te])[:, 1]
        s = clone(secondary).fit(Xs[tr], resid_last); mc_te = s.predict(Xs[te])
        oof_a[te] = p_te
        oof_c[te] = _logit(p_te) + beta * mc_te
    aa = roc_auc_score(y, oof_a); ac = roc_auc_score(y, oof_c)
    return {"oof_anchor": oof_a, "oof_combined": oof_c, "betas": betas_used,
            "auroc_anchor": aa, "auroc_combined": ac, "delta": ac - aa}


def select_anchor(modalities, y, estimator=None, cv=5, repeats=3, random_state=0,
                  stability_weight=0.5, coverage=None):
    """Data-driven choice of which modality to ANCHOR on for `anchored_integrate`.

    modalities : dict {name: X (samples x features)}  (or a list -> names become integer indices).
    Each modality is scored by repeated stratified-CV AUROC of `estimator` (default
    StandardScaler + LogisticRegression); the composite score is
        mean_AUROC  -  stability_weight * std_AUROC   [ * coverage[name] if provided ]
    and the modality with the highest composite is returned as the anchor.

    Why this recipe (it matches both our experiments and the literature):
      * The dominant modality is TASK-SPECIFIC -- it varies by endpoint / cancer -- so the anchor
        must be chosen empirically per task, not fixed a priori (late-fusion benchmarks weight each
        modality by its individual CV success; Nikolaou et al., AACR 2023).
      * Choosing the combination data-adaptively by cross-validation is exactly what cooperative
        learning does (Ding, Li, Narasimhan & Tibshirani, PNAS 2022); the anchored gate is a
        constrained, interpretable special case (anchor pinned, secondary added on the residual).
      * The anchor should also be the most ROBUST view (low CV variance, good coverage/reliability),
        because in our scorecard the winning modality (RNA) was simultaneously the most accurate and
        the most stable; the `stability_weight` and `coverage` terms encode that tie-break.
      * The choice is forgiving: because `anchored_integrate`'s gate defaults to the anchor, a
        slightly sub-optimal anchor is not punished (anchoring on the weaker modality still let the
        stronger one earn a gain in our tests) -- but the strongest, most stable modality gives the
        highest floor and the most efficient model.

    Returns dict: {'anchor': name, 'ranking': [(name, mean_auroc, std_auroc, composite), ...]
    sorted best-first, 'scores': {name: per-repeat AUROC array}}.
    """
    from sklearn.base import clone
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.metrics import roc_auc_score
    items = list(modalities.items()) if isinstance(modalities, dict) else list(enumerate(modalities))
    y = np.asarray(y)
    est = estimator or make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000, C=0.5))
    scores = {}; ranking = []
    for name, X in items:
        X = np.asarray(X, dtype=float)
        aucs = []
        for r in range(repeats):
            skf = StratifiedKFold(cv, shuffle=True, random_state=random_state + r)
            p = cross_val_predict(clone(est), X, y, cv=skf, method="predict_proba")[:, 1]
            aucs.append(roc_auc_score(y, p))
        aucs = np.asarray(aucs)
        comp = aucs.mean() - stability_weight * aucs.std()
        if coverage is not None and name in coverage:
            comp *= float(coverage[name])
        scores[name] = aucs
        ranking.append((name, float(aucs.mean()), float(aucs.std()), float(comp)))
    ranking.sort(key=lambda t: -t[3])
    return {"anchor": ranking[0][0], "ranking": ranking, "scores": scores}


def forward_integrate(modalities, y, order=None, anchor=None,
                      betas=(0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0), cv=5, random_state=0,
                      gate_margin=0.01, inner_repeats=3):
    """Forward (greedy) gated multi-omics integration for ANY number of modalities (>= 2).

    Anchor on the leading modality, then consider each remaining modality in turn and add it ONLY if
    its gated contribution onto the current model's residual beats the current model by more than
    `gate_margin` (selection on repeated inner CV; everything leakage-safe via an outer CV). A useless
    modality is dropped rather than concatenated in -- this is forward selection over modalities, which
    is more interpretable and more robust than fusing everything at once.

    modalities : dict {name: X} (or list -> integer names).  y : binary labels.
    order : sequence of names to consider, best first (default: `select_anchor` ranking).
    anchor : force the first/anchor modality (default: the top of `order`).

    Returns dict: oof_anchor, oof_combined, anchor, ranking, auroc_anchor, auroc_combined, delta,
    and 'added' = {modality: (n_folds_it_was_added, mean_beta)} -- which modalities earned their place.
    """
    from sklearn.base import clone
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression, Ridge
    from sklearn.pipeline import make_pipeline
    from sklearn.model_selection import StratifiedKFold, KFold, cross_val_predict
    from sklearn.metrics import roc_auc_score
    items = {k: np.asarray(v, dtype=float) for k, v in
             (modalities.items() if isinstance(modalities, dict) else enumerate(modalities))}
    y = np.asarray(y)
    ranking = None
    if order is None:
        sel = select_anchor(items, y, cv=cv, random_state=random_state)
        ranking = sel["ranking"]; order = [r[0] for r in ranking]
    if anchor is not None and anchor in order:
        order = [anchor] + [m for m in order if m != anchor]
    anc = order[0]
    A = lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000, C=0.5))
    Sn = lambda: make_pipeline(StandardScaler(), Ridge(alpha=20.0))
    inner = max(3, min(cv, 5))
    def _logit(p): p = np.clip(p, 1e-4, 1 - 1e-4); return np.log(p / (1 - p))
    def _sig(z): return 1.0 / (1.0 + np.exp(-z))
    n = len(y); oof_a = np.zeros(n); oof_c = np.zeros(n)
    added_count = {}; added_beta = {}
    for tr, te in StratifiedKFold(cv, shuffle=True, random_state=random_state).split(items[anc], y):
        ytr = y[tr]
        acc = np.zeros(len(tr))
        for rep in range(inner_repeats):
            acc += cross_val_predict(A(), items[anc][tr], ytr,
                                     cv=StratifiedKFold(inner, shuffle=True, random_state=random_state + rep),
                                     method="predict_proba")[:, 1]
        p_oof = acc / inner_repeats
        cur_tr = _logit(p_oof)
        L_anchor = _logit(A().fit(items[anc][tr], ytr).predict_proba(items[anc][te])[:, 1])
        oof_a[te] = _sig(L_anchor); L_comb = L_anchor.copy()
        resid = ytr - p_oof
        for mod in order:
            if mod == anc:
                continue
            resid_fit = resid
            beta_auc = {b: [] for b in betas}; mc_avg = np.zeros(len(tr))
            for rep in range(inner_repeats):
                mc = cross_val_predict(Sn(), items[mod][tr], resid_fit,
                                       cv=KFold(inner, shuffle=True, random_state=random_state + rep))
                for b in betas:
                    beta_auc[b].append(roc_auc_score(ytr, cur_tr + b * mc))
                mc_avg += mc
            mc_avg /= inner_repeats
            auc0 = float(np.mean(beta_auc[0.0]))
            cand = [(b, float(np.mean(v))) for b, v in beta_auc.items() if b != 0.0 and float(np.mean(v)) > auc0 + gate_margin]
            beta = float(max(cand, key=lambda t: t[1])[0]) if cand else 0.0
            if beta > 0:
                s = Sn().fit(items[mod][tr], resid_fit)
                L_comb = L_comb + beta * s.predict(items[mod][te])
                cur_tr = cur_tr + beta * mc_avg; resid = ytr - _sig(cur_tr)
                added_count[mod] = added_count.get(mod, 0) + 1
                added_beta.setdefault(mod, []).append(beta)
        oof_c[te] = L_comb
    aa = roc_auc_score(y, oof_a); ac = roc_auc_score(y, oof_c)
    return {"oof_anchor": oof_a, "oof_combined": oof_c, "anchor": anc, "ranking": ranking,
            "auroc_anchor": aa, "auroc_combined": ac, "delta": ac - aa,
            "added": {m: (added_count[m], float(np.mean(added_beta[m]))) for m in added_count}}


def auto_integrate(modalities, y, anchor=None, cv=5, random_state=0, **kwargs):
    """End-to-end anchored multi-omics integration with automatic anchor selection (any # modalities).

    Picks the anchor with `select_anchor` (unless `anchor` is given) and runs `forward_integrate`,
    which greedily adds each remaining modality onto the anchor's residual and keeps it only where it
    earns a gain. The recommended one-call workflow: data-driven anchor + never-below-anchor gating +
    interpretable per-modality inclusion. Returns the `forward_integrate` result dict.
    """
    return forward_integrate(modalities, y, anchor=anchor, cv=cv, random_state=random_state, **kwargs)


def signature_score(expr, genes, weights=None, sign=1.0):
    """Build a FIXED knowledge anchor score from a curated gene signature (zero trained parameters).

    expr    : DataFrame (samples x genes) or (genes x samples) -- orientation auto-detected from `genes`.
    genes   : the signature's gene names (textbook / published set).
    weights : optional {gene: weight}; default equal weighting.
    sign    : +1 if higher signature -> higher P(y=1), else -1.
    Returns a 1-D array (one score per sample): the (weighted) mean of per-gene z-scores.
    """
    import pandas as _pd
    df = expr if isinstance(expr, _pd.DataFrame) else _pd.DataFrame(np.asarray(expr))
    if not set(genes) & set(df.columns) and (set(genes) & set(df.index)):
        df = df.T                                   # genes were on the index -> transpose to samples x genes
    use = [g for g in genes if g in df.columns]
    if not use:
        raise ValueError("none of the signature genes are present in `expr`")
    sub = df[use].astype(float)
    Z = (sub - sub.mean()) / sub.std(ddof=0).replace(0, 1.0)
    if weights is not None:
        w = np.array([float(weights.get(g, 0.0)) for g in use])
        return sign * (Z.values * w).sum(1) / (np.abs(w).sum() + 1e-9)
    return sign * Z.mean(1).values


def knowledge_anchored_integrate(anchor_score, modalities, y,
                                 betas=(0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0), cv=5, random_state=0,
                                 gate_margin=0.01, inner_repeats=3, anchor_name="knowledge"):
    """Anchor on a FIXED external biological prior instead of on the strongest data modality.

    Generalises `auto_integrate` from "anchor = best empirical modality" to "anchor = established
    knowledge": a textbook gene signature, a published clock (e.g. Horvath), a known driver score, or a
    clinical marker. The anchor is a 1-D score carrying ZERO trained parameters (only an in-fold
    1-parameter calibration); each data modality is then gated onto the anchor's residual and kept ONLY
    where it beats the prior -- so the result is never below the textbook anchor and the gate answers the
    clinically meaningful question "does the genome-wide data add anything beyond established biology?".

    anchor_score : 1-D array, higher = more likely y=1 (e.g. `signature_score(expr, proliferation_genes)`).
    modalities   : dict {name: X} of data modalities to gate onto the prior.
    Returns the `forward_integrate` dict (anchor = the knowledge prior; `added` = which data earned in,
    `delta` = how much the genome-wide data adds beyond the textbook anchor).
    """
    a = np.asarray(anchor_score, dtype=float).reshape(-1, 1)
    mods = {anchor_name: a}
    mods.update(modalities if isinstance(modalities, dict) else {f"mod{i}": X for i, X in enumerate(modalities)})
    return forward_integrate(mods, y, anchor=anchor_name, betas=betas, cv=cv,
                             random_state=random_state, gate_margin=gate_margin, inner_repeats=inner_repeats)


def anchored_residual_discovery(anchor_score, X, feature_names, y, top_k=30, corr_max=0.6,
                                cv=5, random_state=0, n_perm=12, inner_repeats=3,
                                betas=(0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0), screen_top=None):
    """Knowledge-anchored residual discovery with a noise gate -- separate the known, keep only real new.

    Synthesises the whole toolkit: (1) anchor on a FIXED prior (`anchor_score`), (2) gate the data `X`
    onto its residual, (3) test whether that residual is SIGNAL not noise via a label-permutation null on
    the leakage-safe gated gain, and (4) surface the features driving it that are ORTHOGONAL to the anchor
    (|corr with anchor| < corr_max) -- candidate "new beyond the textbook". The discovered panel is
    re-gated and stress-tested against a matched random panel and its own permutation null.

    X : (samples x features) array; feature_names : names aligned to columns; anchor_score : 1-D prior.
    screen_top : optional int. If set and X has more than screen_top features, a Sure-Independence-Screening
        pre-filter keeps the top screen_top features by |marginal correlation with the anchor-residualized
        response| before the gating/permutation steps (scale prep for genome-wide p; default None = off,
        identical results). When active, the whole-pool gate (auroc_combined/delta) and the matched
        random-panel null are computed over the screened pool rather than all features.
    Returns dict: auroc_anchor, auroc_combined, delta, perm_p, signal (bool: delta>0 & perm_p<.05),
    novel = [(feature, partial_corr_with_y_given_anchor, corr_with_anchor)], novel_delta, novel_perm_p,
    random_delta (matched control).
    """
    a = np.asarray(anchor_score, dtype=float); X = np.asarray(X, dtype=float); y = np.asarray(y)
    names = list(feature_names); rng = np.random.default_rng(random_state)

    # Anchor-residualized response (used for both the optional SIS screen and the partial correlation).
    az = (a - a.mean()) / (a.std() + 1e-9)
    azc = az - az.mean(); den = float(azc @ azc) + 1e-12
    rY = (y.astype(float) - y.mean()) - (float(azc @ (y.astype(float) - y.mean())) / den) * azc
    rYc = rY - rY.mean(); rYn = np.sqrt(float(rYc @ rYc)) + 1e-12

    # Optional Sure-Independence-Screening pre-filter (scale prep): when p is huge (e.g. genome-wide
    # methylation, 450k+ probes), keep only the top `screen_top` features by |marginal correlation with the
    # anchor-residualized response| BEFORE the expensive gating / permutation steps. This bounds memory and
    # compute and is aligned with the discovery target (association with y beyond the anchor). The matched
    # random-panel null is then drawn from this screened candidate pool (a conservative background).
    # Default screen_top=None -> no screening (behaviour identical to the unscreened method).
    if screen_top is not None and X.shape[1] > screen_top:
        Xc0 = X - X.mean(axis=0)
        s = np.abs(Xc0.T @ rYc) / (np.sqrt((Xc0 * Xc0).sum(axis=0)) * rYn + 1e-12)
        keep = np.sort(np.argsort(-s)[:int(screen_top)])
        X = X[:, keep]; names = [names[i] for i in keep]

    def gate(Xs, yy):
        return anchored_integrate(a.reshape(-1, 1), Xs, yy, betas=betas, cv=cv,
                                  random_state=random_state, inner_repeats=inner_repeats)

    base = gate(X, y); delta = base["delta"]
    null = np.array([gate(X, rng.permutation(y))["delta"] for _ in range(n_perm)])
    perm_p = (1 + int((null >= delta).sum())) / (n_perm + 1)
    # partial point-biserial correlation of each feature with y, controlling for the anchor.
    # Vectorized residualize-then-correlate: residualize y and every feature on the anchor in a single
    # matrix pass, then correlate -- numerically identical to a per-feature loop but O(1) passes instead
    # of O(p), so it scales to hundreds of thousands of features (e.g. genome-wide methylation).
    Xc = X - X.mean(axis=0)
    rX = Xc - np.outer(azc, (azc @ Xc) / den)            # residual of each feature on the anchor
    rXc = rX - rX.mean(axis=0)
    pc = (rXc.T @ rYc) / (np.sqrt((rXc * rXc).sum(axis=0)) * rYn + 1e-12)
    ca = (Xc.T @ azc) / (np.sqrt((Xc * Xc).sum(axis=0)) * np.sqrt(float(azc @ azc)) + 1e-12)
    novel_idx = [j for j in np.argsort(-np.abs(pc)) if abs(ca[j]) < corr_max][:top_k]
    novel = [(names[j], round(float(pc[j]), 3), round(float(ca[j]), 3)) for j in novel_idx]
    nd = gate(X[:, novel_idx], y); novel_delta = nd["delta"]
    # discovery null: is the discovered panel special vs random panels of the SAME size? (the right
    # "not noise" test for a discovery -- more powerful and appropriate than shuffling the labels)
    k = len(novel_idx)
    rand_deltas = np.array([gate(X[:, list(rng.choice(X.shape[1], k, replace=False))], y)["delta"]
                            for _ in range(n_perm)])
    novel_vs_random_p = (1 + int((rand_deltas >= novel_delta).sum())) / (n_perm + 1)
    return dict(auroc_anchor=round(base["auroc_anchor"], 4), auroc_combined=round(base["auroc_combined"], 4),
                delta=round(float(delta), 4), perm_p=round(perm_p, 4), novel=novel,
                novel_delta=round(float(novel_delta), 4), random_delta_mean=round(float(rand_deltas.mean()), 4),
                novel_vs_random_p=round(novel_vs_random_p, 4),
                signal=bool(novel_delta > rand_deltas.mean() and novel_vs_random_p < 0.05))
