#!/usr/bin/env python3
"""External validation of auto_integrate on REAL biological subtype labels (TCGA-BRCA).

Runs one endpoint per invocation (keeps each run light) and appends to external_subtype_results.csv.
Endpoints: methylation-defined subtype (methylation_Clusters_nature2012) vs RNA-defined (PAM50).
The question is whether auto_integrate ADAPTS -- methylation should lead the methylation-defined
labels, RNA should lead the expression-defined labels, and the gate should never fall below the leader.

Usage:  python dmoi_external_subtype.py <endpoint>
        endpoints: methC1..methC5 | PAM50Basal | PAM50LumA | PAM50Her2 | PAM50LumB
Env:    FEATURE_DIR (cached rna_mc.tsv/meth_gw.tsv/labels_mc.tsv), RNA_TOPK, METH_TOPK
"""
import os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from omniomics import multiomics as mo

FD = os.environ.get("FEATURE_DIR", "/sessions/sleepy-blissful-allen/mnt/outputs")
RNA_TOPK = int(os.environ.get("RNA_TOPK", "600"))
METH_TOPK = int(os.environ.get("METH_TOPK", "500"))
OUT_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "external_subtype_results.csv")


def load():
    rna = pd.read_csv(f"{FD}/rna_mc.tsv", sep="\t", index_col=0)
    meth = pd.read_csv(f"{FD}/meth_gw.tsv", sep="\t", index_col=0)
    lab = pd.read_csv(f"{FD}/labels_mc.tsv", sep="\t", index_col=0)
    S = [s for s in rna.columns if s in meth.columns and s in lab.index]
    rna = rna[S].loc[rna[S].var(axis=1).sort_values(ascending=False).index[:RNA_TOPK]]
    M = meth[S].T.astype(float); M = M.fillna(M.mean()).fillna(0.0)
    M = M[M.var().sort_values(ascending=False).index[:METH_TOPK]]
    return rna[S].T.values, M.values, lab.reindex(S)


def make_y(lab, ep):
    if ep.startswith("methC"):
        return (lab["methylation_Clusters_nature2012"].values == float(ep[5:])).astype(int)
    if ep.startswith("PAM50"):
        return (lab["PAM50Call_RNAseq"].values == ep[5:]).astype(int)
    raise SystemExit(f"unknown endpoint {ep}")


def main():
    ep = sys.argv[1]
    Xr, Xm, lab = load()
    y = make_y(lab, ep)
    if y.sum() < 20 or (len(y) - y.sum()) < 20:
        print(f"{ep}: skip (imbalanced {y.sum()}/{len(y)})"); return
    t = time.time()
    r = mo.auto_integrate({"RNA": Xr, "methylation": Xm}, y, cv=5, random_state=0, inner_repeats=1)
    kind = "methylation-defined" if ep.startswith("methC") else "expression-defined"
    rk = ", ".join(f"{n}={m:.3f}" for n, m, sd, c in r["ranking"])
    row = dict(endpoint=ep, kind=kind, n=len(y), pos=int(y.sum()), anchor=r["anchor"],
               auroc_anchor=round(r["auroc_anchor"], 4), auroc_combined=round(r["auroc_combined"], 4),
               delta=round(r["delta"], 4), added=";".join(r["added"]) or "none",
               ranking=rk, secs=round(time.time() - t, 1))
    df = pd.DataFrame([row])
    if os.path.exists(OUT_CSV):
        df = pd.concat([pd.read_csv(OUT_CSV), df], ignore_index=True).drop_duplicates("endpoint", keep="last")
    df.to_csv(OUT_CSV, index=False)
    print(f"{ep:12} [{kind}] anchor={r['anchor']:11} | {rk} | combined={r['auroc_combined']:.3f} "
          f"delta={r['delta']:+.3f} added={row['added']}  ({row['secs']}s)")


if __name__ == "__main__":
    main()
