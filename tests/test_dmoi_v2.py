"""Offline unit tests for the DMOI v2 functions in omniomics.multiomics.

No network, no real data — synthetic matrices only. Validates:
  * dmoi_v2_representation column structure (rna_/meth_/int_ per pole) and int == rna*meth
  * reliability weighting actually shifts a pole's methylation toward the high-reliability gene
  * interactions=False / reliability=None degrade gracefully (v1-like main effects)
  * methylation_reliability: correlated probes -> high score, decorrelated -> low, single probe -> NaN
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from omniomics import multiomics as mo


def _toy(seed=0, n=12):
    rng = np.random.default_rng(seed)
    S = [f"s{i}" for i in range(n)]
    genes = ["A", "B", "C"]
    rna = pd.DataFrame(rng.normal(size=(len(genes), n)), index=genes, columns=S)
    meth = pd.DataFrame(rng.normal(size=(len(genes), n)), index=genes, columns=S)
    poles = {"P": ["GS"]}
    gmt = {"GS": ["A", "B", "C"]}
    return rna, meth, poles, gmt, S


def test_v2_columns_and_interaction():
    rna, meth, poles, gmt, S = _toy()
    F = mo.dmoi_v2_representation(rna, meth, poles, gmt, reliability=None, interactions=True)
    assert set(F.columns) == {"rna_P", "meth_P", "int_P"}
    # interaction column is exactly the product of the two main-effect columns
    assert np.allclose(F["int_P"].values, (F["rna_P"] * F["meth_P"]).values)
    # no interaction requested -> only main effects, and never the redundant v1 difference
    F2 = mo.dmoi_v2_representation(rna, meth, poles, gmt, interactions=False)
    assert set(F2.columns) == {"rna_P", "meth_P"}
    assert not any(c.startswith("disagree") or c.startswith("diff") for c in F2.columns)


def test_reliability_weighting_shifts_pole():
    rna, meth, poles, gmt, S = _toy(seed=1)
    # standardized methylation per gene, with meth_sign = -1 baked into the function
    zA = ((-meth.loc["A"]) - (-meth.loc["A"]).mean()) / (-meth.loc["A"]).std(ddof=0)
    rel_on_A = {"A": 1.0, "B": 0.0, "C": 0.0}
    Fw = mo.dmoi_v2_representation(rna, meth, poles, gmt, reliability=rel_on_A)
    Fe = mo.dmoi_v2_representation(rna, meth, poles, gmt, reliability=None)
    # weighting onto A makes the pole-methylation track gene A far more than equal weighting does
    cw = np.corrcoef(Fw["meth_P"].values, zA.values)[0, 1]
    ce = np.corrcoef(Fe["meth_P"].values, zA.values)[0, 1]
    assert cw > 0.99 and cw > ce + 0.1


def test_methylation_reliability_proxy():
    rng = np.random.default_rng(2)
    n = 40
    base = rng.normal(size=n)
    cols = [f"s{i}" for i in range(n)]
    M = pd.DataFrame(
        [
            base + rng.normal(scale=0.01, size=n),   # p1  (gene G1)
            base + rng.normal(scale=0.01, size=n),   # p2  (gene G1) -> ~perfectly correlated
            rng.normal(size=n),                      # p3  (gene G2)
            rng.normal(size=n),                      # p4  (gene G2) -> independent
            rng.normal(size=n),                      # p5  (gene G3, single probe)
        ],
        index=["p1", "p2", "p3", "p4", "p5"], columns=cols,
    )
    g2p = {"G1": ["p1", "p2"], "G2": ["p3", "p4"], "G3": ["p5"]}
    rel = mo.methylation_reliability(M, g2p)
    assert rel["G1"] > 0.8          # consistent probes -> reliable
    assert rel["G2"] < 0.5          # decorrelated probes -> unreliable
    assert rel["G3"] != rel["G3"]   # single probe -> NaN (NaN != NaN)


def test_genelevel_structure_and_shrinkage():
    rna, meth, poles, gmt, S = _toy(seed=3)
    G = mo.dmoi_v2_genelevel(rna, meth, genes=["A", "B", "C"],
                             reliability={"A": 1.0, "B": 0.0, "C": 0.0})
    assert set(G.columns) == {f"{a}_{g}" for a in ("rna", "meth", "int") for g in ("A", "B", "C")}
    assert np.allclose(G["int_A"].values, (G["rna_A"] * G["meth_A"]).values)
    # a zero-reliability gene's methylation feature is shrunk toward 0 relative to a reliable gene
    assert G["meth_B"].abs().mean() < 0.05 * G["meth_A"].abs().mean()


def test_dmoi_regimes_labels():
    df = pd.DataFrame(
        {"occ": [3.0, -3.0, 3.0, 0.1], "meth": [2.0, -2.0, -2.0, 0.1], "expr": [1.0, -1.0, 1.0, 0.1]},
        index=["up", "down", "mixed", "near0"],
    )
    R = mo.dmoi_regimes(df)
    assert {"concordance", "n_up", "regime"}.issubset(R.columns)
    assert R.loc["up", "regime"] == "concordant_up" and R.loc["up", "n_up"] == 3
    assert R.loc["down", "regime"] == "concordant_down" and R.loc["down", "n_up"] == 0
    assert R.loc["mixed", "regime"] == "discordant"


def test_anchored_gate():
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(0); n = 500
    y = (rng.random(n) < 0.5).astype(int)
    anchor = np.clip(0.5 + 0.12 * (2 * y - 1) + rng.normal(0, 0.15, n), 0.01, 0.99)
    a0 = roc_auc_score(y, anchor)
    # an informative secondary -> gate engages (beta>0) and improves on the anchor
    helpful = (2 * y - 1) + rng.normal(0, 0.4, n)
    b1, c1 = mo.anchored_gate(anchor, helpful, y)
    assert b1 > 0 and roc_auc_score(y, c1) > a0
    # any secondary (even pure noise) can never drag the combiner below the anchor (beta=0 fallback)
    noise = rng.normal(0, 1, n)
    b0, c0 = mo.anchored_gate(anchor, noise, y)
    assert roc_auc_score(y, c0) >= a0 - 1e-9


def test_anchored_integrate():
    rng = np.random.default_rng(0); n = 300
    y = (rng.random(n) < 0.5).astype(int)
    sig = (2 * y - 1).astype(float)
    Xa = np.column_stack([0.6 * sig + rng.normal(0, 1, n), rng.normal(0, 1, n)])   # weak anchor signal
    Xs = np.column_stack([2.0 * sig + rng.normal(0, 1, n), rng.normal(0, 1, n)])   # strong orthogonal secondary
    res = mo.anchored_integrate(Xa, Xs, y, cv=5, random_state=0)
    assert res["auroc_combined"] > res["auroc_anchor"] + 0.02   # informative secondary -> real gain
    assert max(res["betas"]) > 0                                # gate engaged in >=1 fold
    # a pure-noise secondary must not meaningfully drag the combiner below the anchor
    res2 = mo.anchored_integrate(Xa, rng.normal(0, 1, (n, 3)), y, cv=5, random_state=0)
    assert res2["auroc_combined"] >= res2["auroc_anchor"] - 0.02


def test_select_anchor():
    rng = np.random.default_rng(0); n = 300
    y = (rng.random(n) < 0.5).astype(int); sig = (2 * y - 1).astype(float)
    strong = np.column_stack([1.5 * sig + rng.normal(0, 1, n), rng.normal(0, 1, n)])
    weak = np.column_stack([0.3 * sig + rng.normal(0, 1, n), rng.normal(0, 1, n)])
    res = mo.select_anchor({"strong": strong, "weak": weak}, y)
    assert res["anchor"] == "strong"
    assert res["ranking"][0][1] > res["ranking"][1][1]   # higher mean AUROC ranks first


def test_auto_integrate():
    rng = np.random.default_rng(1); n = 400
    f1 = rng.normal(size=n); f2 = rng.normal(size=n)
    y = ((f1 + 0.6 * f2) > 0).astype(int)                      # f1 dominant, f2 orthogonal
    A = np.column_stack([f1 + rng.normal(0, 0.3, n), rng.normal(0, 1, n)])   # clean f1 -> anchor
    B = np.column_stack([f2 + rng.normal(0, 0.3, n), rng.normal(0, 1, n)])   # clean f2 -> adds value
    res = mo.auto_integrate({"A": A, "B": B}, y, cv=5, random_state=0)
    assert res["anchor"] == "A"                                # picks the dominant modality
    assert res["auroc_combined"] >= res["auroc_anchor"] - 1e-9 # never below the anchor
    assert res["auroc_combined"] > res["auroc_anchor"]         # orthogonal B earns a real gain
    assert "B" in res["added"]                                 # B was added in >=1 fold


def test_forward_integrate_three():
    rng = np.random.default_rng(2); n = 420
    f1 = rng.normal(size=n); f2 = rng.normal(size=n)
    y = ((f1 + 0.6 * f2) > 0).astype(int)
    A = np.column_stack([f1 + rng.normal(0, 0.3, n), rng.normal(0, 1, n)])   # dominant -> anchor
    B = np.column_stack([f2 + rng.normal(0, 0.3, n), rng.normal(0, 1, n)])   # orthogonal -> added
    C = rng.normal(size=(n, 3))                                              # pure noise -> dropped
    res = mo.forward_integrate({"A": A, "B": B, "C": C}, y, cv=5, random_state=0)
    assert res["anchor"] == "A"
    assert "B" in res["added"] and "C" not in res["added"]     # B earns its place; C is dropped
    assert res["auroc_combined"] > res["auroc_anchor"]


def test_signature_score():
    df = pd.DataFrame({"g1": [1.0, 2, 3, 4, 5], "g2": [5.0, 4, 3, 2, 1], "g3": [0.0, 0, 1, 0, 0]},
                      index=[f"s{i}" for i in range(5)])   # samples x genes
    s = mo.signature_score(df, ["g1", "g2"])
    z = (df[["g1", "g2"]] - df[["g1", "g2"]].mean()) / df[["g1", "g2"]].std(ddof=0)
    assert np.allclose(s, z.mean(1).values)
    assert np.allclose(mo.signature_score(df.T, ["g1", "g2"]), s)        # genes-on-index auto-transposed
    assert np.allclose(mo.signature_score(df, ["g1", "g2"], sign=-1.0), -s)
    assert np.allclose(mo.signature_score(df, ["g1", "g2"], weights={"g1": 1.0, "g2": 0.0}),
                       ((df["g1"] - df["g1"].mean()) / df["g1"].std(ddof=0)).values)


def test_knowledge_anchored_integrate():
    rng = np.random.default_rng(3); n = 400
    f1 = rng.normal(size=n); f2 = rng.normal(size=n)
    y = ((f1 + 0.6 * f2) > 0).astype(int)
    anchor = f1 + rng.normal(0, 0.3, n)                                  # FIXED external prior tracking f1
    data = np.column_stack([f2 + rng.normal(0, 0.3, n), rng.normal(0, 1, n)])   # data carries orthogonal f2
    res = mo.knowledge_anchored_integrate(anchor, {"data": data}, y, cv=5, random_state=0)
    assert res["anchor"] == "knowledge"                                 # the fixed prior is pinned as anchor
    assert res["auroc_combined"] >= res["auroc_anchor"] - 1e-9          # never below the textbook anchor
    assert res["auroc_combined"] > res["auroc_anchor"]                  # orthogonal data earns a gain
    assert "data" in res["added"]


def test_anchored_residual_discovery():
    rng = np.random.default_rng(5); n = 500
    f1 = rng.normal(size=n); f2 = rng.normal(size=n)                    # independent known/new axes
    y = ((f1 + 0.9 * f2) > 0).astype(int)
    anchor = f1 + rng.normal(0, 0.2, n)                                 # fixed prior tracks the KNOWN axis f1
    true = [f2 + rng.normal(0, 0.5, n) for _ in range(5)]              # 5 features carry the NEW axis f2
    noise = [rng.normal(0, 1, n) for _ in range(25)]                   # 25 pure-noise features
    X = np.column_stack(true + noise)
    names = [f"new{i}" for i in range(5)] + [f"noise{i}" for i in range(25)]
    res = mo.anchored_residual_discovery(anchor, X, names, y, top_k=5, n_perm=12, cv=5, random_state=0)
    top = [g for g, _, _ in res["novel"]]
    assert sum(t.startswith("new") for t in top) >= 3                  # recovers the real new-axis features
    assert all(abs(ca) < 0.6 for _, _, ca in res["novel"])            # all anchor-orthogonal
    assert res["novel_delta"] > res["random_delta_mean"]              # discovery beats matched random panels
    assert res["stability"] is None                                   # stability off by default (backward-compat)


def test_residual_discovery_stability():
    """Selection-stability statistic: a complement to the panel-vs-random null that stays informative when
    that null saturates. A planted anchor-orthogonal axis must have stability well above its permuted-label
    null (large positive stability_gain); the permuted null itself is near chance."""
    rng = np.random.default_rng(7); n = 400
    f1 = rng.normal(size=n); f2 = rng.normal(size=n)
    y = ((f1 + 0.9 * f2) > 0).astype(int)
    anchor = f1 + rng.normal(0, 0.2, n)
    true = [f2 + rng.normal(0, 0.5, n) for _ in range(6)]
    noise = [rng.normal(0, 1, n) for _ in range(40)]
    X = np.column_stack(true + noise)
    names = [f"new{i}" for i in range(6)] + [f"noise{i}" for i in range(40)]
    res = mo.anchored_residual_discovery(anchor, X, names, y, top_k=8, n_perm=8, cv=5, random_state=0,
                                         stability_reps=20)
    assert res["stability"] is not None and res["stability_null"] is not None
    assert res["stability_null"] < 0.3                               # permuted-label null near chance
    assert res["stability_gain"] > 0.3                               # real axis recurs well above its null
    assert res["stability"] > res["stability_null"]


def test_marker_correlated_anchor():
    """Data-driven (Venet 2011 meta-PCNA) anchor recipe: the top genes correlated with a canonical marker
    are exactly the ones that co-vary with it -- a reproducible, hand-curation-free anchor definition."""
    import pandas as pd
    rng = np.random.default_rng(0); n, p = 120, 200
    prolif = rng.normal(size=n)
    X = rng.normal(size=(n, p)); X[:, 0] = prolif
    for j in range(1, 11):
        X[:, j] += 1.4 * prolif                       # 10 genes correlated with the marker
    df = pd.DataFrame(X, columns=["PCNA"] + [f"g{j}" for j in range(1, p)])
    genes = mo.marker_correlated_anchor(df, marker="PCNA", top_k=12, exclude_marker=True)
    assert "PCNA" not in genes
    assert sum(g.startswith("g") and int(g[1:]) <= 10 for g in genes) >= 8   # recovers the correlated set


def test_hypothesis_anchor_test():
    """Hypothesis-as-anchor 3-way verdict: a hypothesis orthogonal to the textbook that adds signal is
    SUPPORTED; pure noise is REFUTED; a hypothesis that only echoes the textbook is EXPLAINED_BY_TEXTBOOK."""
    rng = np.random.default_rng(0); n = 500
    a = rng.normal(size=n)                              # textbook axis
    b = rng.normal(size=n)                              # a real, orthogonal novel axis
    y = ((a + b) > 0).astype(int)
    T = a + 0.2 * rng.normal(size=n)                    # textbook anchor tracks a
    supported = b + 0.3 * rng.normal(size=n)            # tracks the orthogonal axis -> should add beyond T
    refuted = rng.normal(size=n)                        # pure noise
    explained = a + 0.3 * rng.normal(size=n)            # echoes the textbook axis -> predicts but redundant
    vs = mo.hypothesis_anchor_test(T, supported, y, cv=4, inner_repeats=1)["verdict"]
    vr = mo.hypothesis_anchor_test(T, refuted, y, cv=4, inner_repeats=1)["verdict"]
    ve = mo.hypothesis_anchor_test(T, explained, y, cv=4, inner_repeats=1)["verdict"]
    assert vs == "SUPPORTED"
    assert vr == "REFUTED"
    assert ve == "EXPLAINED_BY_TEXTBOOK"


if __name__ == "__main__":
    test_v2_columns_and_interaction()
    test_reliability_weighting_shifts_pole()
    test_methylation_reliability_proxy()
    test_genelevel_structure_and_shrinkage()
    test_dmoi_regimes_labels()
    test_anchored_gate()
    test_anchored_integrate()
    test_select_anchor()
    test_auto_integrate()
    print("DMOI v2 unit tests: ALL PASS")
