# Aim 1 — modern reproduction of GSE57575 RNA-seq (nf-core)

Re-derives the paper's RNA-seq from **raw reads** with current best practice, to compare against the
2015 Cufflinks-FPKM results and our n=2 moderated test (`run_golden.py`).

## Dataset (SRA SRP041895 / BioProject PRJNA246728)

8 single-end 101 bp Illumina HiSeq 2000 runs — see `srr_sample_map.csv`:

| SRR | sample | genotype | set |
|-----|--------|----------|-----|
| SRR1283908 | WT_set1 | WT | set1 |
| SRR1283909 | WWD_set1 | WWD | set1 |
| SRR1283910 | R_set1 | R | set1 |
| SRR1283911 | TKO_set1 | TKO | set1 |
| SRR1283912 | WT_set2 | WT | set2 |
| SRR1283913 | WWD_set2 | WWD | set2 |
| SRR1283914 | R_set2 | R | set2 |
| SRR1283915 | TKO_set2 | TKO | set2 |

## Run (on a Linux node with Docker — e.g. the swarm's hrkim-linux, NOT macOS/sandbox)

```bash
bash run_week1.sh        # env check -> nf-core test smoke -> fetchngs (FASTQ) -> nf-core/rnaseq
```

Files:
- `ids.csv` — SRR list for `nf-core/fetchngs` (downloads FASTQ + builds a samplesheet).
- `samplesheet.csv` — hand-written nf-core/rnaseq samplesheet with readable sample names (mirror of
  the fetchngs output; use either).
- `run_week1.sh` — the full Week-1 command sequence (pinned versions; GRCm38/mm10; STAR+Salmon).

Genome note: paper used **mm9**; we use **GRCm38 (mm10, GENCODE M25)** for a modern, well-annotated
reference. Build difference is recorded in the Week-2 "2015 vs 2026" concordance step.

## Outputs → Week 2

`results_rnaseq/star_salmon/salmon.merged.gene_counts.tsv` feeds `../run_modern_de.py`, which runs a
**count-based DESeq2** (Set as a blocking covariate) and writes a concordance table vs the paper /
our moderated test. That closes Aim 1's "2015 vs 2026" deliverable.

## Why this is the modern upgrade
- raw reads → STAR + Salmon (vs TopHat/Bowtie + Cufflinks)
- **counts + DESeq2 negative-binomial** with replicate Set as a covariate (vs FPKM eyeballing / n=2)
- pinned, containerised, reproducible Nextflow run with MultiQC (vs bespoke scripts)
