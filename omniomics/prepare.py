"""Preprocess HM450 methylation into genomic-context matrices used by the BRCA scripts.
Portable (pure Python, no awk/shell). Idempotent: skips outputs that already exist.

Produces, into config.cache_dir():
  pg_promoter.tsv, pg_shore.tsv, pg_enhancer.tsv   (probe -> gene per context)
  promoter_probe_gene.tsv                          (alias of pg_promoter, for run_meth_arm)
  promoter_meth.tsv.gz   (promoter probes,  all samples)
  enh_meth_all.tsv.gz    (enhancer probes,  all samples)
  ctx_meth_lumab.tsv.gz  (all-context probes, LumA/LumB samples)

Run:  python -m omniomics.prepare      (or:  omniomics-prepare-brca)
Requires DMOI_BRCA_DATA to point at dmoi-brca-poc/data (see omniomics.config).
"""
import os, gzip, io, urllib.request
from collections import defaultdict
import bisect
from . import config

HG19_REFGENE = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
HG19_CGI     = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz"

def _dl(url, dest):
    if not os.path.exists(dest) or os.path.getsize(dest) == 0:
        print(f"  downloading {os.path.basename(dest)} ...")
        urllib.request.urlretrieve(url, dest)
    return dest

def _build_probe_context(tcga, cache):
    """Classify HM450 probes by genomic context using hg19 TSS + CpG islands."""
    out = {c: os.path.join(cache, f"pg_{c}.tsv") for c in ("promoter","shore","enhancer")}
    if all(os.path.exists(p) for p in out.values()):
        return out
    rg = _dl(HG19_REFGENE, os.path.join(cache, "refGene_hg19.txt.gz"))
    cg = _dl(HG19_CGI,     os.path.join(cache, "cpgIslandExt_hg19.txt.gz"))
    tss = defaultdict(list)
    with gzip.open(rg, "rt") as f:
        for L in f:
            c = L.rstrip("\n").split("\t")
            chrom, strand, txS, txE, sym = c[2], c[3], int(c[4]), int(c[5]), c[12]
            tss[sym].append((chrom, txS if strand == "+" else txE))
    cgi = defaultdict(list)
    with gzip.open(cg, "rt") as f:
        for L in f:
            c = L.split("\t"); cgi[c[1]].append((int(c[2]), int(c[3])))
    for k in cgi: cgi[k].sort()
    def cgi_dist(chrom, pos):
        arr = cgi.get(chrom)
        if not arr: return 10**9
        starts = [a for a, _ in arr]; i = bisect.bisect_right(starts, pos) - 1; best = 10**9
        for j in (i, i+1):
            if 0 <= j < len(arr):
                s, e = arr[j]; best = min(best, 0 if s <= pos <= e else min(abs(pos-s), abs(pos-e)))
        return best
    fhs = {c: open(out[c], "w") for c in out}
    pm = os.path.join(tcga, "hm450_probemap.tsv")
    with open(pm) as f:
        next(f)
        for L in f:
            pid, gene, chrom, s, e, strand = L.rstrip("\n").split("\t")
            if gene == "." or not s.isdigit(): continue
            pos = int(s); cd = cgi_dist(chrom, pos)
            for g in gene.split(","):
                ts = [t for (gc, t) in tss.get(g, ()) if gc == chrom]
                if not ts: continue
                td = min(abs(pos - t) for t in ts)
                if td <= 1500:               fhs["promoter"].write(f"{pid}\t{g}\n")
                if 0 < cd <= 2000:           fhs["shore"].write(f"{pid}\t{g}\n")
                if td > 5000 and cd > 4000:  fhs["enhancer"].write(f"{pid}\t{g}\n")
    for h in fhs.values(): h.close()
    # alias for run_meth_arm
    import shutil; shutil.copy(out["promoter"], os.path.join(cache, "promoter_probe_gene.tsv"))
    return out

def _probe_set(pg_tsv):
    s = set()
    with open(pg_tsv) as f:
        for L in f: s.add(L.split("\t")[0])
    return s

def _filter_matrix(src_gz, probe_set, out_gz, keep_cols=None, chunk=20000):
    """Stream the HM450 matrix; keep header + rows whose probe is in probe_set; optional column subset."""
    if os.path.exists(out_gz) and os.path.getsize(out_gz) > 0:
        print(f"  {os.path.basename(out_gz)} exists — skip"); return out_gz
    print(f"  filtering -> {os.path.basename(out_gz)} ({len(probe_set)} probes) ...")
    with gzip.open(src_gz, "rt") as fi, gzip.open(out_gz, "wt") as fo:
        header = fi.readline().rstrip("\n").split("\t")
        if keep_cols is not None:
            idx = [0] + [i for i, h in enumerate(header) if i > 0 and h in keep_cols]
            fo.write("\t".join(header[i] for i in idx) + "\n")
        else:
            idx = None; fo.write("\t".join(header) + "\n")
        n = 0
        for line in fi:
            p = line.find("\t")
            if line[:p] in probe_set:
                if idx is None: fo.write(line)
                else:
                    parts = line.rstrip("\n").split("\t"); fo.write("\t".join(parts[i] for i in idx) + "\n")
                n += 1
        print(f"    kept {n} rows")
    return out_gz

def main():
    cache = config.cache_dir(); tcga = config.brca_tcga_dir()
    hm450 = os.path.join(tcga, "HumanMethylation450.gz")
    print(f"[prepare] cache = {cache}")
    print(f"[prepare] HM450 = {hm450}")
    if not os.path.exists(hm450):
        raise FileNotFoundError(f"missing {hm450}; set DMOI_BRCA_DATA correctly")
    pg = _build_probe_context(tcga, cache)
    prom = _probe_set(pg["promoter"]); enh = _probe_set(pg["enhancer"])
    union = prom | _probe_set(pg["shore"]) | enh
    # all-sample matrices
    _filter_matrix(hm450, prom, os.path.join(cache, "promoter_meth.tsv.gz"))
    _filter_matrix(hm450, enh,  os.path.join(cache, "enh_meth_all.tsv.gz"))
    # LumA/LumB-sample union matrix
    import pandas as pd
    coh = pd.read_csv(os.path.join(tcga, "cohort_v2.tsv"), sep="\t")
    lumab = set(coh[coh["group"].isin(["LumA","LumB"])]["sample_id"])
    _filter_matrix(hm450, union, os.path.join(cache, "ctx_meth_lumab.tsv.gz"), keep_cols=lumab)
    print("[prepare] done. Cache ready for the BRCA scripts.")

if __name__ == "__main__":
    main()
