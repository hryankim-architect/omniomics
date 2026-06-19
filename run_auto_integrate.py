#!/usr/bin/env python3
"""Runner for anchored multi-omics integration (omniomics.multiomics.auto_integrate) on TCGA-BRCA.

Self-contained: builds an RNA panel (top-variance genes) and a genome-wide methylation panel from
the raw TCGA matrices (cached to DMOI_FEATURE_DIR), then runs `auto_integrate` on real endpoints and
writes auto_integrate_results.csv. Demonstrates the recommended one-call workflow: data-driven anchor
selection + never-below-anchor forward gating. Mirrors run_dmoi_v2.py. Methods: reports/DMOI_method_assessment.md.

Run:
    python run_auto_integrate.py
Env overrides:
    BRCA_DIR            dir with HiSeqV2.gz, HumanMethylation450.gz, cohort_v2.tsv (default: config.brca_tcga_dir())
    DMOI_FEATURE_DIR    where cached rna1500.tsv / meth_gw.tsv live or will be written (default: script dir)
Outputs:
    auto_integrate_results.csv   endpoint, n, chosen anchor, anchor/combined AUROC, delta, modalities added
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

BRCA = os.environ.get("BRCA_DIR", DEFAULT_BRCA)
FD = os.environ.get("DMOI_FEATURE_DIR", os.path.dirname(os.path.abspath(__file__)))


def _lumab():
    coh = pd.read_csv(os.path.join(BRCA, "cohort_v2.tsv"), sep="\t")
    coh = coh[(coh.group.isin(["LumA", "LumB"])) & (coh.has_rna) & (coh.has_meth)]
    return dict(zip(coh.sample_id, coh.group))


def prepare():
    """Load cached RNA/methylation panels, or regenerate them from the raw TCGA matrices."""
    rna_p = os.path.join(FD, "rna1500.tsv"); meth_p = os.path.join(FD, "meth_gw.tsv")
    if os.path.exists(rna_p) and os.path.exists(meth_p):
        return pd.read_csv(rna_p, sep="\t", index_col=0), pd.read_csv(meth_p, sep="\t", index_col=0)
    assert BRCA and os.path.isdir(BRCA), f"BRCA dir not found: {BRCA!r} (set BRCA_DIR)"
    samps = list(_lumab())
    hdr = pd.read_csv(os.path.join(BRCA, "HiSeqV2.gz"), sep="\t", nrows=0).columns
    cols = [c for c in hdr if c == "sample" or c in set(samps)]
    rna = pd.read_csv(os.path.join(BRCA, "HiSeqV2.gz"), sep="\t", index_col=0, usecols=cols)
    rna = rna.loc[rna.var(axis=1).sort_values(ascending=False).index[:1500]]
    rna.to_csv(rna_p, sep="\t")
    with gzip.open(os.path.join(BRCA, "HumanMethylation450.gz"), "rt") as fh:
        samples = fh.readline().rstrip("\n").split("\t")[1:]; rows = {}
        for i, line in enumerate(fh):
            if i % 160 == 0:
                t = line.rstrip("\n").split("\t"); rows[t[0]] = t[1:]
    meth = pd.DataFrame(rows, index=samples).T.apply(pd.to_numeric, errors="coerce")
    meth.to_csv(meth_p, sep="\t")
    return rna, meth


def main():
    rna, meth = prepare()
    lab = _lumab()
    S = [s for s in rna.columns if s in meth.columns and s in lab]
    Xr = rna[S].T.values
    Mm = meth[S].T.astype(float); Mm = Mm.fillna(Mm.mean()).fillna(0.0)
    rows = []

    # endpoint 1 -- LumA/B (RNA-defined): auto_integrate should anchor on RNA and add nothing
    y = np.array([1 if lab[s] == "LumB" else 0 for s in S])
    r = mo.auto_integrate({"RNA": Xr, "methylation": Mm.values}, y, cv=5, random_state=0, inner_repeats=2)
    rows.append(["LumA_vs_LumB", len(S), r["anchor"], r["auroc_anchor"], r["auroc_combined"], r["delta"],
                 ";".join(r["added"]) or "none"])

    # endpoint 2 -- methylation-defined positive control: a held-out methylation axis RNA cannot see
    cps = list(Mm.columns); h = len(cps) // 2
    yA = (Mm[cps[:h]].mean(1).values > np.median(Mm[cps[:h]].mean(1).values)).astype(int)
    r2 = mo.auto_integrate({"RNA": Xr, "methylation_B": Mm[cps[h:]].values}, yA, cv=5, random_state=0, inner_repeats=2)
    rows.append(["methylation_defined_posctrl", len(S), r2["anchor"], r2["auroc_anchor"], r2["auroc_combined"],
                 r2["delta"], ";".join(r2["added"]) or "none"])

    out = pd.DataFrame(rows, columns=["endpoint", "n", "anchor", "auroc_anchor", "auroc_combined", "delta", "added"])
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "auto_integrate_results.csv")
    out.to_csv(p, index=False)
    print(out.round(4).to_string(index=False))
    print(f"\nwrote {p}")
    print("note: auto_integrate picks the dominant modality per task (RNA for LumA/B, methylation for the")
    print("      methylation-defined endpoint) and never falls below it; see reports/DMOI_method_assessment.md.")


if __name__ == "__main__":
    main()
