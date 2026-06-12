#!/usr/bin/env bash
# Install the nf-core toolchain into a DEDICATED conda env (never the base env — installing
# Nextflow's deps into base downgrades numpy/pandas and breaks the Python analysis env).
# Nextflow is a JVM tool and does not share the Python environment.
set -euo pipefail

ENV=${NF_ENV:-nf}
echo "== creating conda env '$ENV' with Java 17 + Nextflow =="
conda create -y -n "$ENV" -c bioconda -c conda-forge "openjdk=17" nextflow

# Container engine: Docker if usable; otherwise rootless Apptainer (add to the same env).
if command -v docker >/dev/null && docker info >/dev/null 2>&1; then
  echo "Docker available -> PROFILE=docker (default)"
else
  echo "No usable Docker -> adding Apptainer to env '$ENV' ; run with PROFILE=singularity"
  conda install -y -n "$ENV" -c conda-forge apptainer || echo "  (ask admin to install apptainer/singularity if this fails)"
  grep -q NXF_SINGULARITY_CACHEDIR ~/.bashrc 2>/dev/null || \
    echo 'export NXF_SINGULARITY_CACHEDIR=$HOME/.nf_singularity_cache' >> ~/.bashrc
fi

echo ""
echo "Done. Run the pipeline from the '$ENV' env (keeps your Python env untouched):"
echo "    conda run -n $ENV bash run_week1.sh                 # Docker"
echo "    conda run -n $ENV env PROFILE=singularity bash run_week1.sh   # HPC"
