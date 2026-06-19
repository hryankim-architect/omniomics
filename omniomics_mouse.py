#!/usr/bin/env python3
"""Shared mm10 (GRCm38) genomic-context annotation for the GSE57577 RRBS + ChIP downstream steps.

Lazily downloads UCSC refGene + cpgIslandExt for mm10, builds merged interval sets for
  promoter  (TSS +/- 2 kb)
  cgi       (CpG island)
  cgi_shore (within 2 kb of a CGI but outside it)
  distal    (everything else; an enhancer-like proxy — see NOTE)
and classifies single positions (CpGs) or intervals (ChIP peaks) in O(log n).

NOTE on "enhancer": mm10 lacks a single canonical enhancer track here, so distal-intergenic is
used as an enhancer-like proxy and is labelled `distal` honestly. Swap in a real enhancer BED
(e.g. ENCODE cCRE) by passing extra intervals to build_contexts(extra={"enhancer": [...]}).

The download is guarded: importing this module does nothing on the network. Geometry functions
(merge_intervals, in_intervals, interval_hits_label, classify_*) are pure and unit-tested offline.
"""
import bisect, gzip, os, sys, urllib.request

PROMOTER_FLANK = 2000
SHORE_FLANK = 2000
UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database"


def cache_dir():
    here = os.path.dirname(os.path.abspath(__file__))
    d = os.environ.get("OMNIOMICS_MM10_CACHE", os.path.join(here, "data", "mm10_cache"))
    os.makedirs(d, exist_ok=True)
    return d


# ---------- pure geometry (offline-testable) ----------
def merge_intervals(ivals):
    """Merge a list of (start, end) into sorted, non-overlapping intervals."""
    s = sorted((int(a), int(b)) for a, b in ivals if b > a)
    out = []
    for a, b in s:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def _starts(merged):
    return [a for a, _ in merged]


def in_intervals(merged, starts, pos):
    """Is point `pos` inside any merged interval? O(log n)."""
    i = bisect.bisect_right(starts, pos) - 1
    return i >= 0 and pos < merged[i][1]


def interval_overlaps(merged, starts, a, b):
    """Does [a, b) overlap any merged interval? O(log n)."""
    i = bisect.bisect_right(starts, b) - 1
    while i >= 0 and merged[i][1] > a:
        if merged[i][0] < b and merged[i][1] > a:
            return True
        i -= 1
        if i >= 0 and merged[i][1] <= a:
            break
    return False


def subtract(a_merged, b_merged):
    """Interval set difference a \\ b (both merged, same chrom)."""
    out, bi = [], 0
    for a0, a1 in a_merged:
        cur = a0
        while bi < len(b_merged) and b_merged[bi][1] <= a0:
            bi += 1
        j = bi
        while j < len(b_merged) and b_merged[j][0] < a1:
            b0, b1 = b_merged[j]
            if b0 > cur:
                out.append((cur, min(b0, a1)))
            cur = max(cur, b1)
            j += 1
            if cur >= a1:
                break
        if cur < a1:
            out.append((cur, a1))
    return merge_intervals(out)


def _normchr(c):
    """Normalize a chromosome name to UCSC style (the mm10 context tracks are UCSC-named).
    iGenomes/Ensembl GRCm38 emits '1','2','X','MT'; UCSC uses 'chr1','chr2','chrX','chrM'.
    Without this, Ensembl-named cov/peak files match no intervals and everything falls to 'distal'."""
    c = str(c)
    if c.startswith("chr"):
        return c
    if c in ("MT", "M"):
        return "chrM"
    return "chr" + c


class Contexts:
    """Per-chromosome merged interval sets with cached start arrays for fast lookup."""

    def __init__(self, label_to_chrom_ivals):
        # label -> {chrom -> merged intervals}
        self.data = {}
        self.starts = {}
        for label, chrom_ivals in label_to_chrom_ivals.items():
            self.data[label] = {c: merge_intervals(v) for c, v in chrom_ivals.items()}
            self.starts[label] = {c: _starts(m) for c, m in self.data[label].items()}

    def classify_point(self, chrom, pos, priority=("promoter", "cgi", "cgi_shore")):
        chrom = _normchr(chrom)
        for label in priority:
            m = self.data.get(label, {}).get(chrom)
            if m and in_intervals(m, self.starts[label][chrom], pos):
                return label
        return "distal"

    def classify_interval(self, chrom, a, b, priority=("promoter", "cgi", "cgi_shore")):
        chrom = _normchr(chrom)
        for label in priority:
            m = self.data.get(label, {}).get(chrom)
            if m and interval_overlaps(m, self.starts[label][chrom], a, b):
                return label
        return "distal"

    def label_set(self, chrom, a, b):
        """All labels an interval touches (for co-localisation / multi-membership)."""
        chrom = _normchr(chrom)
        hits = set()
        for label, chroms in self.data.items():
            m = chroms.get(chrom)
            if m and interval_overlaps(m, self.starts[label][chrom], a, b):
                hits.add(label)
        return hits or {"distal"}


# ---------- UCSC download + builder (network; guarded) ----------
def _fetch(table):
    dest = os.path.join(cache_dir(), f"{table}.txt.gz")
    if not os.path.exists(dest):
        url = f"{UCSC}/{table}.txt.gz"
        print(f"[mm10] downloading {url}")
        urllib.request.urlretrieve(url, dest)
    return dest


def build_contexts():
    """Download (once) and assemble promoter/cgi/cgi_shore interval sets for mm10."""
    from collections import defaultdict
    prom = defaultdict(list)
    refgene = _fetch("refGene")
    with gzip.open(refgene, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            # refGene: bin,name,chrom,strand,txStart,txEnd,...
            chrom, strand, txStart, txEnd = f[2], f[3], int(f[4]), int(f[5])
            if "_" in chrom:  # skip alt/random contigs
                continue
            tss = txStart if strand == "+" else txEnd
            prom[chrom].append((tss - PROMOTER_FLANK, tss + PROMOTER_FLANK))

    cgi = defaultdict(list)
    cpg = _fetch("cpgIslandExt")
    with gzip.open(cpg, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            # cpgIslandExt: bin,chrom,chromStart,chromEnd,name,...
            chrom, s, e = f[1], int(f[2]), int(f[3])
            if "_" in chrom:
                continue
            cgi[chrom].append((s, e))

    cgi_m = {c: merge_intervals(v) for c, v in cgi.items()}
    shore = {}
    for c, m in cgi_m.items():
        flank = merge_intervals([(s - SHORE_FLANK, e + SHORE_FLANK) for s, e in m])
        shore[c] = subtract(flank, m)

    return Contexts({"promoter": prom, "cgi": cgi_m, "cgi_shore": shore})


if __name__ == "__main__":
    ctx = build_contexts()
    n = {lab: sum(len(v) for v in ch.values()) for lab, ch in ctx.data.items()}
    print("[mm10] context interval counts:", n)
