#!/usr/bin/env python3
"""Is cross-omics DISAGREEMENT signal or noise? Three tests on gene-matched RNA + methylation.

DMOI's premise is that where RNA and methylation *disagree* carries information. The linear difference
z_RNA - z_meth is provably inert (it lies in the span of the two main effects), so the only honest test
uses the nonlinear interaction z_RNA * z_meth. On TCGA-BRCA LumA/B (250 variance-top genes with >=3 CpGs,
n=417), we run:

  1. PAIRING-PERMUTATION NULL (the cleanest test). Fit a nonlinear model on [RNA, meth] (mains) vs
     [RNA, meth, INT] where INT = RNA*meth. The increment = AUROC(full) - AUROC(mains). Then break ONLY
     the cross-omics link: rebuild INT from each patient's RNA paired with a *random* patient's meth,
     keeping the real mains. If the true increment exceeds this permuted-pairing null, the interaction's
     value requires the genuine RNA-meth pairing -> it is real joint signal, not a feature artefact.
  2. SYNERGY / interaction-information. Reduce each modality to a 1-D cross-fitted score; the
     information-theoretic synergy S = I(y; R, M) - I(y; R) - I(y; M) (co-information; S>0 = synergy)
     plus a model synergy = AUROC(nonlinear joint) - AUROC(additive), each with a label-permutation null.
  3. HELD-OUT REPLICATION. Repeated train/test splits: does the interaction increment reproduce
     out-of-sample (not just in CV)?

We also report the LINEAR difference increment as a negative control (expected ~0). Writes
discordance_test_results.csv. Honest reading of prior work: the interaction signal is real but small
(~+0.01 AUROC); these tests quantify whether that small signal survives the noise null.

Run:  BRCA_DIR=/path/to/tcga_brca python reports/dmoi_discordance_tests.py
Env:  BRCA_DIR, FG_DIR (cache of disc_rna_z.tsv/disc_meth_z.tsv/disc_y.tsv), DT_PERMS (default 19),
      DT_SPLITS (default 10).
"""
import os, sys, json, gzip
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import cross_val_predict, StratifiedKFold, train_test_split
from sklearn.metrics import roc_auc_score, mutual_info_score

HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
FG = os.environ.get("FG_DIR", os.getcwd())
BRCA = os.environ.get("BRCA_DIR", "")
PERMS = int(os.environ.get("DT_PERMS", "19")); SPLITS = int(os.environ.get("DT_SPLITS", "10"))
RNG = np.random.default_rng(0)


def _hgb():
    return HistGradientBoostingClassifier(max_iter=120, max_depth=3, learning_rate=0.1, random_state=0)


def _oof(X, y, est):
    p = cross_val_predict(est, X, y, cv=StratifiedKFold(5, shuffle=True, random_state=0),
                          method="predict_proba")[:, 1]
    return roc_auc_score(y, p)


def load_or_build():
    rp, mp, yp = (os.path.join(FG, f) for f in ("disc_rna_z.tsv", "disc_meth_z.tsv", "disc_y.tsv"))
    if all(os.path.exists(p) for p in (rp, mp, yp)):
        R = pd.read_csv(rp, sep="\t", index_col=0); M = pd.read_csv(mp, sep="\t", index_col=0)
        y = pd.read_csv(yp, sep="\t", index_col=0).iloc[:, 0]
        g = [x for x in R.index if x in M.index]; S = list(y.index)
        return R.loc[g, S].T.values, M.loc[g, S].T.values, y.values.astype(int)
    # build gene-matched RNA + gene-level methylation from raw (250 top-variance genes with >=3 CpGs)
    assert BRCA and os.path.isdir(BRCA), f"BRCA dir not found: {BRCA!r} (set BRCA_DIR)"
    coh = pd.read_csv(os.path.join(BRCA, "cohort_v2.tsv"), sep="\t")
    coh = coh[(coh.group.isin(["LumA", "LumB"])) & (coh.has_rna) & (coh.has_meth)]
    lab = dict(zip(coh.sample_id, coh.group))
    with gzip.open(os.path.join(BRCA, "HiSeqV2.gz"), "rt") as fh:
        rcols = set(fh.readline().rstrip("\n").split("\t")[1:])
    S = [s for s in coh.sample_id if s in rcols]
    pm = pd.read_csv(os.path.join(BRCA, "hm450_probemap.tsv"), sep="\t", usecols=["#id", "gene"])
    g2c = {}
    for cg, gs in zip(pm["#id"], pm["gene"]):
        if gs == "." or pd.isna(gs):
            continue
        for g in str(gs).split(","):
            g2c.setdefault(g, []).append(cg)
    ge3 = {g for g, cs in g2c.items() if len(cs) >= 3}
    rna = pd.read_csv(os.path.join(BRCA, "HiSeqV2.gz"), sep="\t", index_col=0, usecols=["sample"] + S)
    cand = [g for g in rna.index if g in ge3]
    panel = list(rna.loc[cand].var(axis=1).sort_values(ascending=False).index[:250])
    Rg = rna.loc[panel, S]
    Rz = Rg.sub(Rg.mean(1), axis=0).div(Rg.std(1).replace(0, 1), axis=0)
    want = set(c for g in panel for c in g2c[g])
    keep = {}
    with gzip.open(os.path.join(BRCA, "HumanMethylation450.gz"), "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")[1:]
        for line in fh:
            cg, _, rest = line.partition("\t")
            if cg in want:
                keep[cg] = rest.rstrip("\n").split("\t")
    cm = pd.DataFrame(keep, index=hdr).T.apply(pd.to_numeric, errors="coerce").reindex(columns=S)
    Mrows = {g: cm.reindex([c for c in g2c[g] if c in cm.index]).mean(0) for g in panel}
    Mg = pd.DataFrame(Mrows).T.reindex(columns=S)
    Mg = Mg.apply(lambda r: r.fillna(r.mean()), axis=1).fillna(0.0)
    Mz = Mg.sub(Mg.mean(1), axis=0).div(Mg.std(1).replace(0, 1), axis=0)
    common = [g for g in Rz.index if g in Mz.index]
    y = np.array([1 if lab[s] == "LumB" else 0 for s in S])
    return Rz.loc[common, S].T.values, Mz.loc[common, S].T.values, y


def test_pairing_null(R, M, y):
    INT = R * M
    main = _oof(np.hstack([R, M]), y, _hgb())
    full = _oof(np.hstack([R, M, INT]), y, _hgb())
    inc = full - main
    nulls = []
    for _ in range(PERMS):
        perm = RNG.permutation(len(y))
        nulls.append(_oof(np.hstack([R, M, R * M[perm]]), y, _hgb()) - main)
    nulls = np.array(nulls)
    p = (1 + int((nulls >= inc).sum())) / (PERMS + 1)
    # linear-difference negative control (logistic, additive)
    lr = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
    lin_main = _oof(np.hstack([R, M]), y, lr); lin_diff = _oof(np.hstack([R, M, R - M]), y, lr)
    return dict(test="pairing_permutation", auroc_main=round(main, 4), auroc_full_int=round(full, 4),
                increment_int=round(inc, 4), null_mean=round(float(nulls.mean()), 4),
                p_value=round(p, 4), linear_diff_increment=round(lin_diff - lin_main, 4))


def test_synergy(R, M, y):
    lr = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
    skf = StratifiedKFold(5, shuffle=True, random_state=0)
    rs = cross_val_predict(lr, R, y, cv=skf, method="predict_proba")[:, 1]
    ms = cross_val_predict(lr, M, y, cv=skf, method="predict_proba")[:, 1]
    Z = np.column_stack([rs, ms])
    add = _oof(Z, y, make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)))
    joint = _oof(Z, y, _hgb())
    # information-theoretic synergy on tercile-discretised scores (bits)
    r3 = np.digitize(rs, np.quantile(rs, [1/3, 2/3])); m3 = np.digitize(ms, np.quantile(ms, [1/3, 2/3]))
    nats2bits = 1.0 / np.log(2)
    def S_of(yy):
        return (mutual_info_score(yy, r3 * 3 + m3) - mutual_info_score(yy, r3) - mutual_info_score(yy, m3)) * nats2bits
    S = S_of(y)
    null = np.array([S_of(RNG.permutation(y)) for _ in range(200)])
    pS = (1 + int((null >= S).sum())) / 201
    return dict(test="synergy", auroc_rna=round(roc_auc_score(y, rs), 4), auroc_meth=round(roc_auc_score(y, ms), 4),
                auroc_additive=round(add, 4), auroc_joint_nonlinear=round(joint, 4),
                model_synergy=round(joint - add, 4), info_synergy_bits=round(float(S), 4),
                info_synergy_p=round(pS, 4))


def test_heldout(R, M, y):
    INT = R * M; incs = []
    for k in range(SPLITS):
        idx = np.arange(len(y))
        tr, te = train_test_split(idx, test_size=0.4, stratify=y, random_state=k)
        m = _hgb().fit(np.hstack([R, M])[tr], y[tr]); f = _hgb().fit(np.hstack([R, M, INT])[tr], y[tr])
        am = roc_auc_score(y[te], m.predict_proba(np.hstack([R, M])[te])[:, 1])
        af = roc_auc_score(y[te], f.predict_proba(np.hstack([R, M, INT])[te])[:, 1])
        incs.append(af - am)
    incs = np.array(incs)
    return dict(test="heldout_replication", splits=SPLITS, mean_increment=round(float(incs.mean()), 4),
                frac_positive=round(float((incs > 0).mean()), 2), sd=round(float(incs.std()), 4))


def main():
    R, M, y = load_or_build()
    print(f"gene-matched LumA/B: n={len(y)} genes={R.shape[1]} LumB={int(y.sum())}")
    rows = [test_pairing_null(R, M, y), test_synergy(R, M, y), test_heldout(R, M, y)]
    for r in rows:
        print("  " + "  ".join(f"{k}={v}" for k, v in r.items()))
    # long-format CSV
    recs = []
    for r in rows:
        for k, v in r.items():
            if k != "test":
                recs.append(dict(test=r["test"], metric=k, value=v))
    out = pd.DataFrame(recs); p = os.path.join(REPO, "discordance_test_results.csv"); out.to_csv(p, index=False)
    print("\nwrote", p)


if __name__ == "__main__":
    main()
