#!/usr/bin/env bash
# Aim 1 — modern reproduction of GSE57577 ChIP-seq from raw reads with nf-core/chipseq.
# Run on a Linux node with Docker (e.g. the swarm's hrkim-linux), from the 'nf' conda env:
#   conda run -n nf bash run_chipseq.sh
#
# Dataset: 10 single-end 51 bp HiSeq-2000 runs (SRR1283896–SRR1283906, no 902), SRP041894 / PRJNA246729.
#   Dnmt3a2 occupancy ChIP: WT/WWD/R x (clone-1 + bulk)  -> the engineered-ADD redistribution signal.
#   H3K4me3 ChIP: WT/WWD/R clone-1.   Single input: TKO_input (used as control for all — see NOTE).
#
# NOTE: the study deposited ONE input (TKO genotype). nf-core needs a control per ChIP, so all ChIPs
#       point at TKO_input. That is biologically imperfect (input is from Dnmt-TKO); an analyst may
#       prefer per-genotype inputs if generating new data. Recorded here for transparency.
set -euo pipefail
cd "$(dirname "$0")"

export NXF_SYNTAX_PARSER=v1     # nf-core legacy config compat on Nextflow 25/26
export NXF_ANSI_LOG=false       # SSH-friendly streaming logs

CHIPSEQ_VER=2.1.0
FETCHNGS_VER=1.12.0
GENOME=GRCm38                   # paper used mm9; GRCm38/mm10 for a modern reference
PROFILE=${PROFILE:-docker}      # override: PROFILE=singularity bash run_chipseq.sh
THREADS=${THREADS:-8}
MEM=${MEM:-32.GB}
READLEN=${READLEN:-50}          # reads are 51 bp; chipseq 2.1.0 only allows 50/75/100/150/200 -> 50

echo "== Step 1: environment check (PROFILE=$PROFILE) =="
command -v nextflow >/dev/null || { echo "Nextflow missing -> run ../install_tools.sh"; exit 1; }
command -v java >/dev/null     || { echo "Java missing -> conda install -c conda-forge 'openjdk=17'"; exit 1; }
case "$PROFILE" in
  docker)      command -v docker >/dev/null || { echo "Docker missing -> use PROFILE=singularity"; exit 1; } ;;
  singularity) command -v singularity >/dev/null || command -v apptainer >/dev/null || { echo "Singularity/Apptainer missing"; exit 1; } ;;
esac
nextflow -version | head -1

echo "== Step 2: smoke test (test profile — MUST pass before the real run) =="
nextflow run nf-core/chipseq -r $CHIPSEQ_VER -profile test,$PROFILE --outdir test_chipseq

echo "== Step 3: fetch raw FASTQ from SRA (nf-core/fetchngs) =="
nextflow run nf-core/fetchngs -r $FETCHNGS_VER -profile $PROFILE --input ids.csv --outdir sra

echo "== Step 4: build the chipseq samplesheet (antibody/control) from the downloaded FASTQ =="
python ../build_samplesheet.py --mode chipseq --srdir sra --map srr_sample_map.csv --out samplesheet.csv

echo "== Step 5: nf-core/chipseq (BWA + MACS2) =="
# chipseq 2.1.0 still uses the legacy --max_cpus/--max_memory (check_max). Pin them to the ACTUAL
# node size (auto-detected) — the THREADS=8 default exceeds a 4-CPU node and would fail.
# Override to coexist with another running pipeline, e.g.  CPUS=2 MEM_GB=10 bash run_chipseq.sh
CPUS=${CPUS:-$(nproc)}
MEM_GB=${MEM_GB:-$(free -g 2>/dev/null | awk '/^Mem:/{print $2}')}; [ -z "$MEM_GB" ] && MEM_GB=8; [ "$MEM_GB" -lt 4 ] && MEM_GB=4
echo "  node: ${CPUS} CPUs, ${MEM_GB} GB RAM"
nextflow run nf-core/chipseq -r $CHIPSEQ_VER -profile $PROFILE \
    --input samplesheet.csv \
    --outdir results_chipseq \
    --genome $GENOME \
    --read_length $READLEN \
    --max_cpus $CPUS --max_memory "${MEM_GB}.GB" \
    -resume

echo "== Step 6: key outputs =="
echo "  peaks:     results_chipseq/bwa/mergedLibrary/macs2/narrowPeak/  (Dnmt3a2 occupancy, H3K4me3)"
echo "  bigwigs:   results_chipseq/bwa/mergedLibrary/bigwig/"
echo "  QC:        results_chipseq/multiqc/"
echo "Next: compare Dnmt3a2 occupancy redistribution (WT vs WWD vs R) at promoters/enhancers vs the paper."
