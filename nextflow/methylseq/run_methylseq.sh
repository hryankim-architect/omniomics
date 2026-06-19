#!/usr/bin/env bash
# Aim 1 — modern reproduction of GSE57577 RRBS methylation from raw reads with nf-core/methylseq.
# Run on a Linux node with Docker (e.g. the swarm's hrkim-linux), from the 'nf' conda env:
#   conda run -n nf bash run_methylseq.sh
#
# Dataset: 9 single-end 51 bp HiSeq-2000 RRBS runs (SRR1283918–SRR1283926), SRP041896 / PRJNA246727.
#   WT/WWD/R x2 reps + C(parental) x2 + TKO.  C = unmodified parental ESC control; TKO = Dnmt-TKO.
set -euo pipefail
cd "$(dirname "$0")"

export NXF_SYNTAX_PARSER=v1     # nf-core legacy config compat on Nextflow 25/26
export NXF_ANSI_LOG=false       # SSH-friendly streaming logs

METHYLSEQ_VER=4.2.0
FETCHNGS_VER=1.12.0
GENOME=GRCm38                   # paper used mm9; GRCm38/mm10 for a modern reference (note in concordance)
PROFILE=${PROFILE:-docker}      # override: PROFILE=singularity bash run_methylseq.sh
THREADS=${THREADS:-8}
MEM=${MEM:-32.GB}

echo "== Step 1: environment check (PROFILE=$PROFILE) =="
command -v nextflow >/dev/null || { echo "Nextflow missing -> run ../install_tools.sh"; exit 1; }
command -v java >/dev/null     || { echo "Java missing -> conda install -c conda-forge 'openjdk=17'"; exit 1; }
case "$PROFILE" in
  docker)      command -v docker >/dev/null || { echo "Docker missing -> use PROFILE=singularity"; exit 1; } ;;
  singularity) command -v singularity >/dev/null || command -v apptainer >/dev/null || { echo "Singularity/Apptainer missing"; exit 1; } ;;
esac
nextflow -version | head -1

echo "== Step 2: smoke test (test profile — MUST pass before the real run) =="
nextflow run nf-core/methylseq -r $METHYLSEQ_VER -profile test,$PROFILE --outdir test_methylseq

echo "== Step 3: fetch raw FASTQ from SRA (nf-core/fetchngs) =="
nextflow run nf-core/fetchngs -r $FETCHNGS_VER -profile $PROFILE --input ids.csv --outdir sra

echo "== Step 4: build the methylseq samplesheet from the downloaded FASTQ =="
python ../build_samplesheet.py --mode methylseq --srdir sra --map srr_sample_map.csv --out samplesheet.csv

echo "== Step 5: nf-core/methylseq (Bismark, RRBS mode) =="
# methylseq 4.x DROPPED --max_cpus/--max_memory; cap resources via process.resourceLimits, auto-
# detected from this node (otherwise a process requests more CPUs than exist and the run fails).
# CPUS/MEM_GB auto-detect from the node, but can be overridden to coexist with another running
# pipeline, e.g.  CPUS=2 MEM_GB=10 bash run_methylseq.sh   (leaves headroom for a concurrent run).
CPUS=${CPUS:-$(nproc)}
MEM_GB=${MEM_GB:-$(free -g 2>/dev/null | awk '/^Mem:/{print $2}')}; [ -z "$MEM_GB" ] && MEM_GB=8; [ "$MEM_GB" -lt 4 ] && MEM_GB=4
cat > resources.config <<CFG
process {
    resourceLimits = [ cpus: ${CPUS}, memory: '${MEM_GB}.GB', time: '240.h' ]
}
CFG
echo "  node: ${CPUS} CPUs, ${MEM_GB} GB RAM -> wrote resources.config"
# --rrbs disables dedup and applies RRBS-appropriate trimming (MspI digest). Bismark/bowtie2 on the
# mouse genome fits well under 24 GB.
nextflow run nf-core/methylseq -r $METHYLSEQ_VER -profile $PROFILE \
    --input samplesheet.csv \
    --outdir results_methylseq \
    --genome $GENOME \
    --rrbs \
    -c resources.config \
    -resume

echo "== Step 6: key outputs =="
echo "  per-cytosine methylation: results_methylseq/bismark/methylation_calls/methylation_coverage/*.cov.gz"
echo "  reports:                  results_methylseq/multiqc/"
echo "Next: build genotype x context (promoter/CGI-shore/enhancer) methylation matrix vs the paper's RRBS."
