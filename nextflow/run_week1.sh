#!/usr/bin/env bash
# Aim 1 / Week 1 — modern reproduction of GSE57575 RNA-seq from raw reads with nf-core.
# Run on a Linux node with Docker (e.g. the swarm's hrkim-linux), NOT on macOS/sandbox.
#
# Steps: env check -> nf-core test profile (smoke) -> fetch raw FASTQ from SRA -> nf-core/rnaseq.
# Dataset: 8 single-end 101 bp HiSeq-2000 runs (SRR1283908–SRR1283915 = WT/WWD/R/TKO × set1/2).
set -euo pipefail
cd "$(dirname "$0")"

# Nextflow 25/26 default to a new strict config parser that rejects nf-core 3.22.2's
# legacy `if (...) { process {...} }` config syntax. Fall back to the v1 parser.
export NXF_SYNTAX_PARSER=v1
# Stream the log line-by-line (the default in-place ANSI log looks "frozen" over SSH).
export NXF_ANSI_LOG=false

RNASEQ_VER=3.22.2
FETCHNGS_VER=1.12.0
GENOME=GRCm38            # modern mouse build (mm10; GENCODE M25). Paper used mm9 — note in concordance table.
PROFILE=${PROFILE:-docker}   # override e.g.  PROFILE=singularity bash run_week1.sh
THREADS=${THREADS:-8}
MEM=${MEM:-32.GB}

echo "== Step 1: environment check (PROFILE=$PROFILE) =="
command -v nextflow >/dev/null || { echo "Nextflow missing -> run: bash install_tools.sh   (or: conda install -c bioconda -c conda-forge nextflow)"; exit 1; }
command -v java >/dev/null || { echo "Java missing -> conda install -c conda-forge 'openjdk=17'"; exit 1; }
case "$PROFILE" in
  docker)      command -v docker >/dev/null || { echo "Docker missing -> install Docker, or use PROFILE=singularity"; exit 1; } ;;
  singularity) command -v singularity >/dev/null || command -v apptainer >/dev/null || { echo "Singularity/Apptainer missing -> conda install -c conda-forge apptainer"; exit 1; } ;;
esac
nextflow -version | head -1

echo "== Step 2: smoke test (test profile, tiny data — MUST pass before real run) =="
nextflow run nf-core/rnaseq -r $RNASEQ_VER -profile test,$PROFILE --outdir test_rnaseq
ls test_rnaseq/multiqc/*/multiqc_report.html && echo "SMOKE OK"

echo "== Step 3: fetch raw FASTQ from SRA (nf-core/fetchngs) =="
# ids.csv = one SRR per line. Produces ./fastq/ + an auto samplesheet.
nextflow run nf-core/fetchngs -r $FETCHNGS_VER -profile $PROFILE \
    --input ids.csv --outdir sra --nf_core_pipeline rnaseq
# fetchngs writes FASTQ under sra/fastq/ and a ready samplesheet at sra/samplesheet/samplesheet.csv
# (our hand-written samplesheet.csv mirrors it with readable sample names.)

echo "== Step 4: nf-core/rnaseq (count-based) =="
# Modern nf-core (rnaseq 3.14+) DROPPED --max_cpus/--max_memory. Cap resources via a
# process.resourceLimits config instead, auto-detected from this node. Without this, TrimGalore
# requests 12 CPUs and fails on a small node ("req: 12; avail: 4").
CPUS=$(nproc)
MEM_GB=$(free -g 2>/dev/null | awk '/^Mem:/{print $2}')
[ -z "$MEM_GB" ] && MEM_GB=8
[ "$MEM_GB" -lt 4 ] && MEM_GB=4
cat > resources.config <<CFG
process {
    resourceLimits = [ cpus: ${CPUS}, memory: '${MEM_GB}.GB', time: '240.h' ]
}
CFG
echo "  node: ${CPUS} CPUs, ${MEM_GB} GB RAM -> wrote resources.config"

# STAR alignment of the mouse genome loads the index into ~30 GB RAM. On a small node, default to
# Salmon selective-alignment quantification only (low RAM, gives gene counts) — the count matrix is
# equivalent in spirit to STAR+Salmon for our DESeq2 concordance. Force STAR with ALIGNER=star_salmon
# only if the node has >= 32 GB.
ALIGNER=${ALIGNER:-auto}
if [ "$ALIGNER" = "auto" ]; then
  if [ "$MEM_GB" -ge 32 ]; then ALIGNER=star_salmon; else ALIGNER=salmon_only; fi
fi
if [ "$ALIGNER" = "salmon_only" ]; then
  QUANT="--pseudo_aligner salmon --skip_alignment"
  COUNTS_DIR=salmon
  echo "  quant: Salmon only (--skip_alignment) — low-RAM path"
else
  QUANT="--aligner star_salmon"
  COUNTS_DIR=star_salmon
  echo "  quant: STAR + Salmon"
fi

nextflow run nf-core/rnaseq -r $RNASEQ_VER -profile $PROFILE \
    --input sra/samplesheet/samplesheet.csv \
    --outdir results_rnaseq \
    --genome $GENOME \
    $QUANT \
    -c resources.config \
    -resume

echo "== Step 5: key outputs =="
echo "  counts: results_rnaseq/${COUNTS_DIR}/salmon.merged.gene_counts.tsv"
echo "  tpm:    results_rnaseq/${COUNTS_DIR}/salmon.merged.gene_tpm.tsv"
echo "Next (Week 2): DESeq2 on the count matrix vs the paper's FPKM / our n=2 moderated test."
