#!/usr/bin/env python3
"""Runner for DMOI v2 (reliability-weighted, interaction-aware multi-omics fusion).

Wires the library API (omniomics.multiomics.dmoi_v2_representation + methylation_reliability,
and the v1 dmoi_representation for comparison) into a real benchmark on TCGA-BRCA LumA-vs-LumB.

Self-contained: it regenerates gene-level methylation AND per-gene reliability straight from the
raw HumanMethylation450 matrix + probe map, so it does NOT depend on any derived cache that may
have been cleaned up. Methods + results are written up in reports/DMOI_method_assessment.md.

Run:
    python run_dmoi_v2.py
Env overrides:
    BRCA_DIR             dir with HiSeqV2.gz, HumanMethylation450.gz, cohort_v2.tsv, hm450_probemap.tsv
                         (defaults to omniomics.config.brca_tcga_dir())
    DMOI_METH_SUBSET     optional pre-filtered methylation TSV (probe x sample) to skip the 450K stream
Outputs:
    dmoi_v2_auroc.csv    benchmark table (orig DMOI vs v2 variants vs RNA baseline)
"""
import os, sys, gzip
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from omniomics import multiomics as mo
try:
    from omniomics import config
    DEFAULT_BRCA = config.brca_tcga_dir()
except Exception:
    DEFAULT_BRCA = ""

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from scipy import stats

BRCA = os.environ.get("BRCA_DIR") or DEFAULT_BRCA
MARK = {"ER":     ["ESR1", "GATA3", "FOXA1", "XBP1", "PGR"],
        "PROLIF": ["MKI67", "CCNB1", "CCNE1", "BUB1", "AURKA", "MYBL2", "CDK1", "E2F1", "FOXM1"]}
POLES = {"ER": ["ER"], "PROLIF": ["PROLIF"]}            # pole -> gene-set name(s)
GMT = {k: v for k, v in MARK.items()}                  # gene-set name -> genes
POLEGENES = [g for v in MARK.values() for g in v]


def stream_rows(path, want):
    """{key in `want`} -> list[str] values, plus the sample (column) header, from a gzipped matrix."""
    with gzip.open(path, "rt") as fh:
        samples = fh.readline().rstrip("\n").split("\t")[1:]
        rows = {}
        for line in fh:
            i = line.find("\t")
            if line[:i] in want:
                rows[line[:i]] = line.rstrip("\n").split("\t")[1:]
    return samples, rows


def load_methylation_probes(target_probes):
    """probes x samples DataFrame for `target_probes`; honours DMOI_METH_SUBSET if set."""
    sub = os.environ.get("DMOI_METH_SUBSET")
    if sub and os.path.exists(sub):
        df = pd.read_csv(sub, sep="\t", index_col=0)
        return df.loc[[p for p in df.index if p in target_probes]]
    samples, rows = stream_rows(os.path.join(BRCA, "HumanMethylation450.gz"), target_probes)
    return pd.DataFrame({p: pd.to_numeric(pd.Series(v, index=samples), errors="coerce")
                         for p, v in rows.items()}).T


def cv_auroc(make_model, X, y, gbt=False, rep=12):
    out = []
    for r in range(rep):
        a = []
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=r).split(X, y):
            if gbt:
                clf = make_model().fit(X[tr], y[tr]); pr = clf.predict_proba(X[te])[:, 1]
            else:
                sc = StandardScaler().fit(X[tr])
                clf = make_model().fit(sc.transform(X[tr]), y[tr]); pr = clf.predict_proba(sc.transform(X[te]))[:, 1]
            a.append(roc_auc_score(y[te], pr))
        out.append(np.mean(a))
    return np.array(out)


def main():
    assert BRCA and os.path.isdir(BRCA), f"BRCA dir not found: {BRCA!r} (set BRCA_DIR)"
    coh = pd.read_csv(os.path.join(BRCA, "cohort_v2.tsv"), sep="\t")
    coh = coh[(coh.group.isin(["LumA", "LumB"])) & (coh.has_rna) & (coh.has_meth)]
    lab = dict(zip(coh.sample_id, coh.group))

    # RNA (gene x sample) for pole genes
    rsamp, rrows = stream_rows(os.path.join(BRCA, "HiSeqV2.gz"), set(POLEGENES))
    rna = pd.DataFrame({g: pd.to_numeric(pd.Series(v, index=rsamp), errors="coerce") for g, v in rrows.items()}).T

    # probe -> gene, then methylation probes -> gene-level + reliability (LIBRARY call)
    g2p = {g: [] for g in POLEGENES}
    with open(os.path.join(BRCA, "hm450_probemap.tsv")) as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            if len(f) >= 2 and f[1] != ".":
                for g in f[1].split(","):
                    if g in g2p:
                        g2p[g].append(f[0])
    target = set(p for ps in g2p.values() for p in ps)
    mprobe = load_methylation_probes(target)
    rel = mo.methylation_reliability(mprobe, g2p)               # <-- library API
    med = np.nanmedian([v for v in rel.values() if v == v])
    rel = {g: (rel[g] if rel.get(g) == rel.get(g) else med) for g in g2p}
    methg = {g: mprobe.loc[[p for p in g2p[g] if p in mprobe.index]].mean(axis=0)
             for g in POLEGENES if any(p in mprobe.index for p in g2p[g])}
    meth = pd.DataFrame(methg).T                                # gene x sample

    # matched samples / labels
    S = [s for s in rna.columns if s in lab and s in meth.columns]
    y = np.array([1 if lab[s] == "LumB" else 0 for s in S])
    rna, meth = rna[S], meth[S]

    # representations via the package
    D1 = mo.dmoi_representation(rna, meth, POLES, GMT)          # v1 (rna/meth/disagree)
    D2 = mo.dmoi_v2_representation(rna, meth, POLES, GMT, reliability=rel)        # v2 (rel-weighted + int)
    D2u = mo.dmoi_v2_representation(rna, meth, POLES, GMT, reliability=None)      # v2 ablation (no reliability)
    Dg = mo.dmoi_v2_genelevel(rna, meth, genes=POLEGENES, reliability=rel)       # v2 gene-resolution (item 3)
    pk = list(POLES)
    feats = {
        "RNA poles only (LR)":            (D1[[f"rna_{p}" for p in pk]].values, False),
        "orig DMOI v1 (LR)":              (D1.values, False),
        "DMOI v2 no-reliability (GBT)":   (D2u.values, True),
        "DMOI v2 reliability+int (GBT)":  (D2.values, True),
        "DMOI v2 gene-level (GBT)":       (Dg.values, True),
    }
    lr = lambda: LogisticRegression(max_iter=5000)
    gb = lambda: GradientBoostingClassifier(n_estimators=100, max_depth=2, random_state=0)
    R = {k: cv_auroc(gb if g else lr, X, y, gbt=g) for k, (X, g) in feats.items()}
    base = R["orig DMOI v1 (LR)"]

    print(f"[run_dmoi_v2] BRCA LumA/B  n={len(S)} (LumB={int(y.sum())})  12x5 CV")
    print(f"  mean methylation reliability (inter-probe corr): {np.nanmean([v for v in rel.values()]):.3f}")
    rows = []
    for k, (X, g) in feats.items():
        d = R[k].mean() - base.mean()
        p = 1.0 if k == "orig DMOI v1 (LR)" else stats.wilcoxon(R[k], base).pvalue
        print(f"  {k:32s} AUROC={R[k].mean():.4f}  dVsV1={d:+.4f}  p={p:.1e}")
        rows.append({"model": k, "AUROC_mean": R[k].mean(), "AUROC_std": R[k].std(),
                     "delta_vs_v1": d, "wilcoxon_p_vs_v1": p})

    # permuted-label null (roadmap item 5b): with shuffled labels the v2 model must collapse to ~0.5,
    # proving the reported gains are not an artefact of the CV/feature-construction pipeline.
    Xbest = feats["DMOI v2 no-reliability (GBT)"][0]
    rngp = np.random.default_rng(0)
    null = float(np.mean([cv_auroc(gb, Xbest, rngp.permutation(y), gbt=True, rep=1).mean()
                          for _ in range(20)]))
    print(f"  {'DMOI v2 (permuted-label null)':32s} AUROC={null:.4f}   (must be ~0.5)")
    rows.append({"model": "DMOI v2 (permuted-label null)", "AUROC_mean": null, "AUROC_std": 0.0,
                 "delta_vs_v1": null - base.mean(), "wilcoxon_p_vs_v1": float("nan")})

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dmoi_v2_auroc.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"  wrote {out}")
    print("  note: reference plain ~1500-gene RNA on this task ~= 0.94 (see dmoi_enhancer_auroc.csv);")
    print("        v2's value is robustness + interpretability, not beating high-dim RNA (see report 6).")


if __name__ == "__main__":
    main()
