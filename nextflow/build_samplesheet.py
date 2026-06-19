#!/usr/bin/env python3
"""Build an nf-core samplesheet for the methylseq (RRBS) or chipseq arm of the GSE57577
reproduction, joining our readable srr_sample_map.csv to the FASTQ that nf-core/fetchngs
actually downloaded (whatever it named them).

Usage:
    python build_samplesheet.py --mode methylseq --srdir methylseq/sra \
        --map methylseq/srr_sample_map.csv --out methylseq/samplesheet.csv
    python build_samplesheet.py --mode chipseq   --srdir chipseq/sra \
        --map chipseq/srr_sample_map.csv   --out chipseq/samplesheet.csv

It matches each map row to a FASTQ by looking for the SRR run accession as a substring of the
fastq paths in <srdir>/samplesheet/samplesheet.csv (the fetchngs output) — and falls back to a
glob of <srdir>/fastq/<SRR>*.fastq.gz. This makes the join robust to fetchngs naming (SRR / SRX /
composite). Emits the exact column layout each pipeline version expects.
"""
import argparse, csv, glob, os, sys


def fetchngs_fastq_index(srdir):
    """run_accession-substring -> (fastq_1, fastq_2) from the fetchngs output samplesheet."""
    idx = []
    sheet = os.path.join(srdir, "samplesheet", "samplesheet.csv")
    if os.path.exists(sheet):
        with open(sheet) as fh:
            for row in csv.DictReader(fh):
                f1 = row.get("fastq_1") or row.get("fastq_1".upper()) or ""
                f2 = row.get("fastq_2") or ""
                blob = " ".join(str(v) for v in row.values())
                idx.append((blob, f1, f2))
    return idx


def resolve_fastq(srr, idx, srdir):
    for blob, f1, f2 in idx:
        if srr in blob:
            return f1, f2
    hits = sorted(glob.glob(os.path.join(srdir, "fastq", f"{srr}*.fastq.gz")))
    if hits:
        f2 = hits[1] if len(hits) > 1 else ""
        return os.path.abspath(hits[0]), (os.path.abspath(f2) if f2 else "")
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["methylseq", "chipseq"])
    ap.add_argument("--srdir", required=True, help="fetchngs output dir (has fastq/ and samplesheet/)")
    ap.add_argument("--map", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    idx = fetchngs_fastq_index(a.srdir)
    rows = list(csv.DictReader(open(a.map)))
    missing = []
    out_rows = []

    for r in rows:
        srr = r["run"]
        f1, f2 = resolve_fastq(srr, idx, a.srdir)
        if not f1:
            missing.append(srr)
            continue
        if a.mode == "methylseq":
            out_rows.append({"sample": r["sample"], "fastq_1": f1, "fastq_2": f2 or ""})
        else:  # chipseq
            out_rows.append({"sample": r["sample"], "fastq_1": f1, "fastq_2": f2 or "",
                             "replicate": r.get("replicate", ""),  # may be blank -> auto-numbered below
                             "antibody": r.get("antibody", ""), "control": r.get("control", ""),
                             "control_replicate": r.get("control_replicate", "")})

    if a.mode == "chipseq":
        # nf-core/chipseq requires a `replicate` column (header is
        #   sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate ).
        # If the map didn't supply one, auto-number replicates per sample in input order, so a
        # sample name that appears N times becomes REP1..REPN (e.g. each Dnmt3a2 ChIP has two
        # runs = clone-1 + bulk -> REP1, REP2). Rows that share sample AND replicate would be
        # merged as technical replicates, which is not what we want here.
        seen = {}
        for row in out_rows:
            if not row.get("replicate"):
                seen[row["sample"]] = seen.get(row["sample"], 0) + 1
                row["replicate"] = str(seen[row["sample"]])

    if missing:
        sys.exit(f"[build-samplesheet] could not locate FASTQ for: {missing}\n"
                 f"  Did fetchngs finish into {a.srdir}/ ?  (looked in samplesheet/ and fastq/)")

    cols = (["sample", "fastq_1", "fastq_2"] if a.mode == "methylseq"
            else ["sample", "fastq_1", "fastq_2", "replicate", "antibody", "control", "control_replicate"])
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(out_rows)
    print(f"[build-samplesheet] wrote {len(out_rows)} samples -> {a.out}")


if __name__ == "__main__":
    main()
