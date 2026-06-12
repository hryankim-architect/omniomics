#!/usr/bin/env python3
"""Aim 1 / Week 2 — count-based DESeq2 reanalysis of the nf-core/rnaseq output, with the replicate
Set as a blocking covariate, and a '2015 vs 2026' concordance vs the paper / our n=2 moderated test.

Input : nextflow/results_rnaseq/star_salmon/salmon.merged.gene_counts.tsv  (produced by run_week1.sh)
Run   : pip install pydeseq2 ; python run_modern_de.py
Output: modern_de_<contrast>.csv + modern_de_concordance.csv
"""
import os, sys, glob
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
COUNTS = os.environ.get("RNASEQ_COUNTS",
    os.path.join(ROOT, "nextflow", "results_rnaseq", "star_salmon", "salmon.merged.gene_counts.tsv"))
MAP = os.path.join(ROOT, "nextflow", "srr_sample_map.csv")

if not os.path.exists(COUNTS):
    sys.exit(f"[modern-de] counts not found: {COUNTS}\n"
             f"  Run nextflow/run_week1.sh on a Docker/Linux node first (Aim 1 Week 1).")
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
except ImportError:
    sys.exit("[modern-de] pip install pydeseq2")

# ---- load Salmon gene counts -> genes x samples (rounded ints) ----
df = pd.read_csv(COUNTS, sep="\t")
gene_col = "gene_name" if "gene_name" in df.columns else df.columns[0]
samp_cols = [c for c in df.columns if c not in ("gene_id", "gene_name")]
counts = df.groupby(gene_col)[samp_cols].sum().round().astype(int).T   # samples x genes

# Resolve each count-matrix column to (genotype, set) regardless of how the pipeline named
# samples. nf-core/fetchngs + rnaseq may label columns by run (SRR), experiment (SRX), GSM, or
# our readable name (WT_set1) — build a lookup over ALL of them and match exactly or by substring
# (handles composite tags like "SRX540720_SRR1283909").
smap = pd.read_csv(MAP)
lut = {}
for _, r in smap.iterrows():
    for col in ("run", "srx", "gsm", "sample"):
        if col in smap.columns and pd.notna(r[col]):
            lut[str(r[col])] = (r["genotype"], r["replicate_set"])

def _resolve(colname):
    c = str(colname)
    if c in lut:
        return lut[c]
    for k, v in lut.items():        # substring match for composite tags
        if k in c:
            return v
    return (None, None)

geno, sett = {}, {}
for c in counts.index:
    geno[c], sett[c] = _resolve(c)
meta = pd.DataFrame({"genotype": pd.Series(geno), "set": pd.Series(sett)}).dropna()

if meta.shape[0] < 8:
    unresolved = [c for c in geno if geno[c] is None]
    sys.exit(f"[modern-de] only resolved {meta.shape[0]}/8 columns to the genotype map.\n"
             f"  count-matrix columns: {list(counts.index)}\n"
             f"  unresolved: {unresolved}\n"
             f"  -> paste the salmon.merged.gene_counts.tsv header and I'll fix the SRR/SRX/GSM map.")

counts = counts.loc[meta.index]
counts = counts.loc[:, counts.sum(axis=0) > 0]
print(f"[modern-de] {counts.shape[0]} samples x {counts.shape[1]} genes; design ~ set + genotype")

# ---- DESeq2 with Set as covariate (the modern upgrade over n=2 FPKM) ----
dds = DeseqDataSet(counts=counts, metadata=meta, design_factors=["set", "genotype"], ref_level=["genotype", "WT"])
dds.deseq2()

GOLDEN = {"WWD": 1888, "R": 9, "TKO": 3}   # our n=2 paired moderated-test counts (FDR<0.05)
TARGETS = ["Gata4", "Dab2", "Lama1", "Col4a1", "Col4a2", "Enc1"]
rows = []
for g in ["WWD", "R", "TKO"]:
    st = DeseqStats(dds, contrast=["genotype", g, "WT"]); st.summary()
    res = st.results_df.dropna(subset=["padj"])
    res.to_csv(os.path.join(ROOT, f"modern_de_{g}_vs_WT.csv"))
    n_sig = int((res["padj"] < 0.05).sum())
    tgt_sig = int(res.reindex(TARGETS).dropna().query("padj < 0.05").shape[0]) if g == "WWD" else None
    rows.append({"contrast": f"{g}_vs_WT", "DESeq2_padj<0.05": n_sig,
                 "our_n2_moderated": GOLDEN[g],
                 "named_targets_sig(/6)": tgt_sig})
con = pd.DataFrame(rows)
con.to_csv(os.path.join(ROOT, "modern_de_concordance.csv"), index=False)
print("\n[2015 vs 2026] concordance:")
print(con.to_string(index=False))
print("\nFull tables: modern_de_*_vs_WT.csv")
