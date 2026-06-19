#!/usr/bin/env bash
# Aim 1 pipeline status — methylseq + chipseq progress at a glance.
#
# Usage:
#   bash check_aim1.sh            # one-shot snapshot
#   watch -n 60 bash check_aim1.sh   # live, refresh every 60s  (Ctrl-C to stop)
#
# Override paths if your layout differs:
#   NF=~/omniomics/nextflow SSD=/mnt/nfssd bash check_aim1.sh
set -uo pipefail
NF="${NF:-$HOME/omniomics/nextflow}"
SSD="${SSD:-/mnt/nfssd}"

hr(){ printf '%.0s-' {1..64}; echo; }

echo "Aim 1 status — $(date '+%Y-%m-%d %H:%M:%S')"
hr

# --- is a pipeline actively running? ---
if pgrep -f 'nf-core/(methylseq|chipseq|fetchngs)' >/dev/null 2>&1; then
  now=$(pgrep -af 'nf-core/(methylseq|chipseq|fetchngs)' | grep -oE 'nf-core/(methylseq|chipseq|fetchngs)' | head -1)
  echo "RUN : RUNNING  (${now})"
elif tmux has-session -t aim1 2>/dev/null; then
  echo "RUN : tmux 'aim1' alive, no nextflow process (between steps, or just finished)"
else
  echo "RUN : not running"
fi

# --- disk (the thing that kept biting us) ---
echo
echo "DISK:"
df -h / "$SSD" 2>/dev/null | sed 1d | awk '{printf "  %-12s size %-5s used %-5s free %-5s (%s)\n",$NF,$2,$3,$4,$5}'

# --- per-arm progress from each pipeline's .nextflow.log ---
for arm in methylseq chipseq; do
  log="$NF/$arm/.nextflow.log"
  echo
  echo "[$arm]"
  if [ ! -f "$log" ]; then echo "  (not started yet)"; continue; fi
  sub=$(grep -c 'Submitted process'  "$log" 2>/dev/null || echo 0)
  comp=$(grep -c 'Task completed'    "$log" 2>/dev/null || echo 0)
  last=$(grep 'Submitted process >'  "$log" 2>/dev/null | tail -1 | sed 's/.*Submitted process > //')
  echo "  tasks: submitted=$sub  completed=$comp"
  echo "  last : ${last:-—}"
  done_ok=$(grep -E 'Pipeline completed successfully' "$log" 2>/dev/null | tail -1)
  [ -n "$done_ok" ] && echo "  >> pipeline completed successfully"
  err=$(grep -E 'ERROR ~|No space left|Pipeline failed|terminated with an error' "$log" 2>/dev/null | tail -1)
  [ -n "$err" ] && echo "  !! ${err}"
done

# --- completion artifacts (the real deliverables) ---
echo
echo "OUTPUTS:"
mcov=$(find "$NF/methylseq/results_methylseq" -name '*.cov.gz' 2>/dev/null | wc -l)
cpk=$(find  "$NF/chipseq/results_chipseq"     -name '*.narrowPeak' 2>/dev/null | wc -l)
printf "  methylseq cov.gz   : %s / 9 %s\n" "$mcov" "$([ "${mcov:-0}" -ge 9 ] && echo '   METHYLSEQ DONE')"
printf "  chipseq narrowPeak : %s %s\n"     "$cpk"  "$([ "${cpk:-0}" -gt 0 ] && echo 'peaks present')"

# --- both arms done? ---
if grep -q 'ALL_DONE' <(tmux capture-pane -t aim1 -p 2>/dev/null); then
  echo
  echo ">>> ALL_DONE marker seen — methylseq + chipseq chain finished."
fi
hr
