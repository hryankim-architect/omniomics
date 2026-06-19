#!/usr/bin/env bash
# Aim 1 FAST path — stop the full-depth methylseq run, downsample the already-downloaded RRBS FASTQ
# to N reads/sample, re-run methylseq on the subsample (genome index reused via -resume), then chipseq.
# Work stays on the ext4 SSD ($SSD).  Downsampling = first-N reads (spatially random on the flowcell,
# so unbiased for genomic methylation), which is plenty for genotype x context MEAN methylation.
#
# Run in tmux (use a NEW session name 'aim1f' — step 0 kills the old 'aim1' run):
#   tmux new-session -d -s aim1f 'bash ~/omniomics/run_aim1_fast.sh; echo ALL_DONE; exec bash'
# Faster/slower: N_READS=5000000 bash run_aim1_fast.sh   (5M ~= a few hours; 15M ~= overnight)
set -u
NF="$HOME/omniomics/nextflow"
SSD="${SSD:-/mnt/nfssd}"
N_READS="${N_READS:-15000000}"
MS="$NF/methylseq"

echo "== 0) stop current full-depth run + nf-core containers =="
tmux kill-session -t aim1 2>/dev/null || true
pkill -9 -f nextflow   2>/dev/null || true
pkill -9 -f 'conda run' 2>/dev/null || true
cids=$(docker ps -q --filter "name=nxf-"); [ -n "$cids" ] && docker kill $cids || true
sleep 3
if pgrep -f 'nf-core/methylseq' >/dev/null; then echo "  !! still running — abort"; exit 1; fi
echo "  stopped."

cd "$MS"
[ -f samplesheet.csv ] || { echo "  !! no samplesheet.csv in $MS"; exit 1; }
echo "== current samplesheet.csv =="; cat samplesheet.csv

echo "== 1) downsample to ${N_READS} reads/sample -> fastq_ds/ =="
mkdir -p "$MS/fastq_ds"
hdr=$(head -1 samplesheet.csv)
i_s=$(echo "$hdr" | tr ',' '\n' | grep -nx 'sample'  | cut -d: -f1)
i_f=$(echo "$hdr" | tr ',' '\n' | grep -nx 'fastq_1' | cut -d: -f1)
[ -n "${i_s:-}" ] && [ -n "${i_f:-}" ] || { echo "  !! can't find sample/fastq_1 columns"; exit 1; }
echo "sample,fastq_1,fastq_2" > samplesheet_ds.csv
while IFS=, read -r sample f1; do
  [ -z "$sample" ] && continue
  src=$(readlink -f "$f1" 2>/dev/null || echo "$f1")
  out="$MS/fastq_ds/${sample}.ds.fastq.gz"
  [ -r "$src" ] || { echo "  !! source not readable: $src"; exit 1; }
  echo "  >> $sample : $(basename "$src")"
  zcat "$src" | head -n $((N_READS*4)) | gzip -1 > "$out"
  [ -s "$out" ] || { echo "  !! downsample produced empty file: $sample"; exit 1; }
  echo "${sample},${out}," >> samplesheet_ds.csv
done < <(awk -F, -v s="$i_s" -v f="$i_f" 'NR>1 && $s!="" {print $s","$f}' samplesheet.csv)
echo "== downsampled samplesheet =="; cat samplesheet_ds.csv

echo "== 2) re-run methylseq on the subsample (work on $SSD) =="
export NXF_WORK="$SSD" NXF_SYNTAX_PARSER=v1 NXF_ANSI_LOG=false
conda run -n nf nextflow run nf-core/methylseq -r 4.2.0 -profile docker \
  --input samplesheet_ds.csv --outdir results_methylseq \
  --genome GRCm38 --rrbs -c resources.config -resume

echo "== 3) chipseq (full depth) =="
conda run -n nf bash "$NF/chipseq/run_chipseq.sh"

echo "== fast path finished =="
