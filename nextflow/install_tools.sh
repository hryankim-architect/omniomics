#!/usr/bin/env bash
# One-shot install of the nf-core toolchain on a Linux node with conda/anaconda (e.g. hrkim-linux).
# Installs Java 17 + Nextflow, and a rootless container engine (Apptainer) for HPC-friendly runs.
set -euo pipefail

echo "== installing Java 17 + Nextflow (+ Apptainer) via conda =="
conda install -y -c bioconda -c conda-forge "openjdk=17" nextflow

# Container engine: Docker if you have it + daemon + docker group; otherwise rootless Apptainer.
if command -v docker >/dev/null && docker info >/dev/null 2>&1; then
  echo "Docker available -> use:  bash run_week1.sh"
else
  echo "No usable Docker -> installing Apptainer (rootless) ; run with PROFILE=singularity"
  conda install -y -c conda-forge apptainer || echo "  (if this fails, ask the admin to install apptainer/singularity)"
  # nf-core reuses the 'singularity' profile for Apptainer; cache images in a roomy dir:
  echo 'export NXF_SINGULARITY_CACHEDIR=$HOME/.nf_singularity_cache' >> ~/.bashrc
fi

echo ""
nextflow -version | head -3 || true
echo "Done. Next:  cd nextflow && PROFILE=${PROFILE:-singularity} bash run_week1.sh"
