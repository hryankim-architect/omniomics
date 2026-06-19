#!/usr/bin/env python3
"""Aim 1 (ChIP-seq) — downstream redistribution analysis for the nf-core/chipseq output.

Quantifies where the engineered Dnmt3a2 binds across genotypes (WT -> WWD -> R) by genomic context
(promoter / CGI / CGI-shore / distal), and how much Dnmt3a2 occupancy co-localises with H3K4me3 —
the modern read-out of the paper's finding that the H3K4me3-recognising ADD (WWD) redistributes
Dnmt3a2 toward active, H3K4me3-marked chromatin.

Input  : nextflow/chipseq/results_chipseq/**/<sample>*_peaks.narrowPeak   (run_chipseq.sh)
Map    : nextflow/chipseq/srr_sample_map.csv   (sample -> genotype, antibody)
Run    : python run_chip_redistribution.py
Output : chip_context_distribution.csv, chip_h3k4me3_colocalization.csv

narrowPeak columns: chrom,start,end,name,score,strand,signalValue,pValue,qValue,peakOffset.
Skip-safe: exits cleanly if no peak files are present yet.
"""
import csv, glob, os, sys
from collections import defaultdict
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get("CHIPSEQ_RESULTS", os.path.join(ROOT, "nextflow", "chipseq", "results_chipseq"))
MAP = os.path.join(ROOT, "nextflow", "chipseq", "srr_sample_map.csv")
CONTEXTS = ["promoter", "cgi", "cgi_shore", "distal"]


def find_peaks():
    return sorted(set(glob.glob(os.path.join(RES, "**", "*.narrowPeak"), recursive=True)))


def read_peaks(path):
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                out.append((f[0], int(f[1]), int(f[2])))
            except ValueError:
                continue
    return out


def sample_of(path, samples):
    base = os.path.basename(path)
    hits = [s for s in samples if s in base]
    return max(hits, key=len) if hits else None


def merged_by_chrom(peaks):
    import omniomics_mouse as mm
    by = defaultdict(list)
    for c, a, b in peaks:
        by[c].append((a, b))
    m = {c: mm.merge_intervals(v) for c, v in by.items()}
    st = {c: [x for x, _ in v] for c, v in m.items()}
    return m, st


def main():
    files = find_peaks()
    if not files:
        sys.exit(f"[chip] no narrowPeak files under {RES}\n"
                 f"  Run nextflow/chipseq/run_chipseq.sh on a Docker/Linux node first.")
    rows = list(csv.DictReader(open(MAP)))
    meta = {r["sample"]: (r["genotype"], r.get("antibody", "")) for r in rows}
    samples = [s for s in meta if meta[s][1]]  # ChIP samples (skip input)

    import omniomics_mouse as mm
    ctx = mm.build_contexts()

    peaks_by_sample = {}
    dist_rows = []
    for path in files:
        samp = sample_of(path, samples)
        if not samp:
            continue
        geno, ab = meta[samp]
        pk = read_peaks(path)
        # A genotype's Dnmt3a2 ChIP is split across REP1/REP2 files that share one sample name —
        # accumulate (don't overwrite) so co-localisation sees the union of both replicates.
        peaks_by_sample.setdefault(samp, []).extend(pk)
        counts = defaultdict(int)
        for c, a, b in pk:
            counts[ctx.classify_interval(c, a, b)] += 1
        tot = max(sum(counts.values()), 1)
        stem = os.path.basename(path).replace("_peaks.narrowPeak", "")  # per-replicate library id
        row = {"library": stem, "sample": samp, "genotype": geno, "antibody": ab, "n_peaks": len(pk)}
        for c in CONTEXTS:
            row["pct_" + c] = round(100.0 * counts[c] / tot, 2)
        dist_rows.append(row)

    if not dist_rows:
        sys.exit("[chip] peak files found but none matched the sample map names.")
    dist = pd.DataFrame(dist_rows).set_index("library")
    dist.to_csv(os.path.join(ROOT, "chip_context_distribution.csv"))

    # H3K4me3 co-localisation: fraction of each genotype's Dnmt3a2 peaks overlapping that
    # genotype's H3K4me3 peaks.
    h3 = {meta[s][0]: peaks_by_sample[s] for s in peaks_by_sample if meta[s][1] == "H3K4me3"}
    coloc_rows = []
    for s, pk in peaks_by_sample.items():
        geno, ab = meta[s]
        if ab != "Dnmt3a2" or geno not in h3:
            continue
        m, st = merged_by_chrom(h3[geno])
        ov = sum(1 for c, a, b in pk if c in m and mm.interval_overlaps(m[c], st[c], a, b))
        coloc_rows.append({"genotype": geno, "dnmt3a2_sample": s, "n_dnmt3a2_peaks": len(pk),
                           "pct_overlapping_H3K4me3": round(100.0 * ov / max(len(pk), 1), 2)})
    coloc = pd.DataFrame(coloc_rows)
    if len(coloc):
        coloc.to_csv(os.path.join(ROOT, "chip_h3k4me3_colocalization.csv"), index=False)

    pd.set_option("display.width", 160)
    print(f"[chip] {len(dist)} ChIP samples")
    print("\n[peak distribution by genomic context (%)]")
    print(dist.to_string())
    if len(coloc):
        print("\n[Dnmt3a2 x H3K4me3 co-localisation]  (expect WWD > R, WT — ADD redistribution)")
        print(coloc.to_string(index=False))
    print("\nWrote chip_context_distribution.csv" + (", chip_h3k4me3_colocalization.csv" if len(coloc) else ""))


if __name__ == "__main__":
    main()
