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

## Install the toolchain (Linux node, once)

nf-core needs **Nextflow + Java 17 + a container engine**. Install into a **dedicated conda env**
(NOT base — installing Nextflow's deps into base downgrades numpy/pandas and breaks the Python env):

```bash
cd nextflow
bash install_tools.sh        # creates conda env 'nf' with Java 17 + Nextflow (+ Apptainer if no Docker)
```

## Run (on the Linux node, NOT macOS/sandbox)

Run the pipeline **from the `nf` env**; keep your Python/omniomics env separate:

```bash
conda run -n nf bash run_week1.sh                              # Docker (default)
conda run -n nf env PROFILE=singularity bash run_week1.sh      # HPC / no Docker
# env check -> nf-core test smoke -> fetchngs (FASTQ) -> nf-core/rnaseq
```

Then run the Python concordance from your omniomics env:
```bash
python ../run_modern_de.py
```

### Troubleshooting
- **"Nextflow missing"** → `bash install_tools.sh` (or `conda install -c bioconda -c conda-forge nextflow`).
- **Java errors** → Nextflow wants Java 11–21; `conda install -c conda-forge 'openjdk=17'`.
- **No Docker on a shared node** → use `PROFILE=singularity`; set `export NXF_SINGULARITY_CACHEDIR=$HOME/.nf_singularity_cache` (roomy disk) so images cache, not re-pull.
- **Memory/CPU** → match the node: `THREADS=16 MEM='64.GB' PROFILE=singularity bash run_week1.sh`.
- **Resume after a failure** → re-run the same command; Nextflow `-resume` continues from the last checkpoint.

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

---

# The other two assays (same study, separate SubSeries)

GSE57577 is a SuperSeries. Each assay was deposited as its own SRA study — all confirmed from the
NCBI SRA RunInfo (single-end 51 bp HiSeq 2000):

| assay | SubSeries SRA | BioProject | runs | dir |
|-------|---------------|------------|------|-----|
| RNA-seq | SRP041895 | PRJNA246728 | SRR1283908–915 (8) | `./` |
| ChIP-seq | SRP041894 | PRJNA246729 | SRR1283896–906, no 902 (10) | `chipseq/` |
| RRBS methylation | SRP041896 | PRJNA246727 | SRR1283918–926 (9) | `methylseq/` |

## RRBS → `methylseq/` (nf-core/methylseq 4.2.0, Bismark, `--rrbs`)

9 runs: WT/WWD/R ×2 reps, C (parental control) ×2, TKO. Run from the `nf` env on a Docker/Linux node:

```bash
conda run -n nf bash methylseq/run_methylseq.sh        # smoke → fetchngs → build sheet → methylseq
```

Outputs per-cytosine coverage (`results_methylseq/bismark/.../*.cov.gz`) → feed a genotype × context
(promoter / CGI-shore / enhancer) methylation matrix to compare against the paper's RRBS and the
HM450-derived contexts used in the DMOI work.

## ChIP-seq → `chipseq/` (nf-core/chipseq 2.1.0, BWA + MACS2)

10 runs: **Dnmt3a2 occupancy** ChIP (WT/WWD/R × clone-1 + bulk — the engineered-ADD redistribution
signal) and **H3K4me3** ChIP (WT/WWD/R clone-1), with a single **TKO_input** control.

```bash
conda run -n nf bash chipseq/run_chipseq.sh            # smoke → fetchngs → build sheet → chipseq
```

> NOTE — the study deposited only **one input** (TKO genotype). nf-core/chipseq needs a control per
> ChIP, so every ChIP points at `TKO_input`. This is biologically imperfect (the input is from
> Dnmt-TKO cells); it is recorded transparently here, and an analyst generating new data would use
> per-genotype inputs. `clone-1` and `bulk` are entered as replicates 1 and 2 of each condition so
> MACS2 consensus peaks can be called.

Outputs narrowPeak + bigWig (`results_chipseq/bwa/mergedLibrary/macs2/`, `.../bigwig/`) → compare
Dnmt3a2 occupancy redistribution (WT→WWD→R) at promoters/enhancers against the paper's ADD-targeting figure.

## Shared helper — `build_samplesheet.py`

Both arms download FASTQ with `nf-core/fetchngs`, then `build_samplesheet.py` joins our readable
`srr_sample_map.csv` to whatever fetchngs named the files (matches on the SRR accession as a
substring; falls back to a `fastq/<SRR>*` glob). It emits the exact column layout each pipeline
version expects (methylseq: `sample,fastq_1,fastq_2`; chipseq adds `antibody,control,control_replicate`)
and **fails loudly**, listing any SRR whose FASTQ it could not locate.

Genome note (all arms): paper used **mm9**; we use **GRCm38 (mm10)** for a modern reference; the
build difference is recorded in each arm's concordance step.
