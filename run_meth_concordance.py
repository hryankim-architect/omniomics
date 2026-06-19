#!/usr/bin/env python3
"""Aim 1 (RRBS) — downstream concordance for the nf-core/methylseq output.

Turns Bismark per-cytosine coverage into a genotype x genomic-context (promoter / CGI / CGI-shore /
distal) methylation table, and the WWD-vs-parental(C) delta — the modern read-out of the paper's
RRBS conclusion that the engineered ADD (WWD) reshapes the methylation landscape.

Input  : nextflow/methylseq/results_methylseq/**/<sample>*.cov.gz   (run_methylseq.sh)
Map    : nextflow/methylseq/srr_sample_map.csv   (sample -> genotype)
Run    : python run_meth_concordance.py
Output : meth_context_by_genotype.csv, meth_context_delta.csv

Bismark .cov.gz columns (no header): chrom, start, end, methylation%, count_meth, count_unmeth.
Skip-safe: exits cleanly if the pipeline output is not present yet.
"""
import csv, glob, gzip, os, sys
from collections import defaultdict
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get("METHYLSEQ_RESULTS", os.path.join(ROOT, "nextflow", "methylseq", "results_methylseq"))
MAP = os.path.join(ROOT, "nextflow", "methylseq", "srr_sample_map.csv")
MIN_COV = int(os.environ.get("RRBS_MIN_COV", "5"))
CONTEXTS = ["promoter", "cgi", "cgi_shore", "distal"]


def find_cov_files():
    files = glob.glob(os.path.join(RES, "**", "*.cov.gz"), recursive=True)
    # bismark also emits *.bismark.cov.gz; keep all coverage files
    return sorted(set(files))


def sample_of(path, samples):
    base = os.path.basename(path)
    hits = [s for s in samples if s in base]
    return max(hits, key=len) if hits else None


def parse_cov_line(line, min_cov=MIN_COV):
    """Bismark .cov line -> (chrom, pos, methylation%) or None if malformed/under-covered."""
    f = line.rstrip("\n").split("\t")
    if len(f) < 6:
        return None
    try:
        chrom, start = f[0], int(f[1])
        meth, unmeth = int(f[4]), int(f[5])
    except ValueError:
        return None
    if meth + unmeth < min_cov:
        return None
    return chrom, start, 100.0 * meth / (meth + unmeth)


def main():
    cov = find_cov_files()
    if not cov:
        sys.exit(f"[meth] no Bismark coverage under {RES}\n"
                 f"  Run nextflow/methylseq/run_methylseq.sh on a Docker/Linux node first.")
    smap = {r["sample"]: r["genotype"] for r in csv.DictReader(open(MAP))}
    samples = list(smap)

    import omniomics_mouse as mm
    ctx = mm.build_contexts()  # downloads mm10 refGene + cpgIslandExt once

    # per (genotype, context): accumulate sum of site methylation% and site count
    msum = defaultdict(float)
    mcnt = defaultdict(int)
    used = []
    for path in cov:
        samp = sample_of(path, samples)
        if not samp:
            continue
        geno = smap[samp]
        used.append((samp, geno, os.path.basename(path)))
        op = gzip.open(path, "rt")
        for line in op:
            site = parse_cov_line(line)
            if site is None:
                continue
            chrom, start, pct = site
            c = ctx.classify_point(chrom, start)
            msum[(geno, c)] += pct
            mcnt[(geno, c)] += 1
        op.close()

    if not used:
        sys.exit("[meth] coverage files found but none matched the sample map names.")
    genos = sorted({g for _, g in [(s, g) for s, g, _ in used]})
    rows = []
    for g in genos:
        row = {"genotype": g}
        for c in CONTEXTS:
            n = mcnt[(g, c)]
            row[c] = (msum[(g, c)] / n) if n else float("nan")
            row[c + "_n"] = n
        rows.append(row)
    by = pd.DataFrame(rows).set_index("genotype")
    by.to_csv(os.path.join(ROOT, "meth_context_by_genotype.csv"))

    # WWD-vs-parental(C) and other-vs-C deltas (mean methylation change per context)
    base = "C" if "C" in by.index else ("WT" if "WT" in by.index else by.index[0])
    delta = by[CONTEXTS].subtract(by.loc[base, CONTEXTS], axis=1).drop(index=base)
    delta.to_csv(os.path.join(ROOT, "meth_context_delta.csv"))

    pd.set_option("display.width", 160)
    print(f"[meth] {len(used)} coverage files, genotypes={genos}, baseline={base}")
    print("\n[mean % methylation by genotype x context]")
    print(by[CONTEXTS].round(2).to_string())
    print(f"\n[delta vs {base}]  (positive = gain of methylation)")
    print(delta.round(2).to_string())
    print("\nWrote meth_context_by_genotype.csv, meth_context_delta.csv")


if __name__ == "__main__":
    main()
