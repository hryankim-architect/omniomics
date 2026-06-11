#!/usr/bin/env bash
# Install omniomics into the Agentic Bioinformatics Swarm on this host (e.g. chi-mac-p).
# Run from the omniomics-prototype repo root:  bash swarm/install.sh
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
DEST="${OMNIOMICS_HOME:-$HOME/omniomics}"
SWARM="${SWARM_HOME:-$HOME/swarm}"

echo "[1/4] copy engine -> $DEST"
mkdir -p "$DEST"
rsync -a --exclude data/ --exclude '__pycache__' --exclude '*.egg-info' "$REPO"/ "$DEST"/

echo "[2/4] editable install"
( cd "$DEST" && pip install -e . )

echo "[3/4] wire Claude Code layer into $SWARM/.claude"
mkdir -p "$SWARM/.claude/commands" "$SWARM/.claude/agents"
cp "$DEST/swarm/commands/golden.md"        "$SWARM/.claude/commands/"
cp "$DEST/swarm/agents/omniomics-runner.md" "$SWARM/.claude/agents/"
echo "    -> merge $DEST/swarm/settings.snippet.json into $SWARM/.claude/settings.json (permissions.allow)"

echo "[4/4] smoke test: golden GSE57577 reproduction"
( cd "$DEST" && python run_golden.py )

cat <<EOF

✅ omniomics installed into the swarm.
   Next (inside 'claude' on the swarm):
     > @omniomics-runner add a dataset to the cohort and harmonize
     > /golden                 # regression gate (blocks merge on FAIL)
     > @swarm-reviewer review the last diff
   For the BRCA pipeline:  export DMOI_BRCA_DATA=~/Downloads/AI/dmoi-brca-poc/data && omniomics-prepare-brca
EOF
