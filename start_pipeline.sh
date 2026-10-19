#!/bin/bash
#SBATCH --job-name=xraft_scrna
#SBATCH --output=slurm/%x.%j.out
#SBATCH --error=slurm/%x.%j.err
#SBATCH --time=144:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=long
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kristian.gurashi@ndcls.ox.ac.uk

# ── Guard: Slurm opens slurm/<job>.out BEFORE this script runs. If slurm/ is
# missing the job is cancelled at 0 s with exit 0:53 and no log at all.
# Running this file directly (./start_pipeline.sh) creates the directory and
# submits the job for you; running sbatch yourself requires slurm/ to exist.
if [ -z "${SLURM_JOB_ID:-}" ]; then
    cd "$(dirname "${BASH_SOURCE[0]}")"
    mkdir -p slurm slurm/sps run
    exec sbatch "${BASH_SOURCE[0]}" "$@"
fi

# ── Environment modules: a batch job is a non-interactive shell, so `module`
# may not be defined (it is only inherited when sbatch is run from a login
# shell). Initialise it explicitly rather than relying on that.
if ! command -v module >/dev/null 2>&1; then
    for f in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh \
             /usr/share/lmod/lmod/init/bash /usr/share/Modules/init/bash \
             /etc/profile; do
        if [ -r "$f" ]; then source "$f" 2>/dev/null && command -v module >/dev/null 2>&1 && break; fi
    done
fi
if ! command -v module >/dev/null 2>&1; then
    echo "[ERROR] 'module' command not available in the job environment" >&2
    exit 1
fi

module load nextflow/25.04.7

# Submit from the pipeline directory after `mkdir -p slurm` (sbatch needs the
# --output directory to exist). tidy_pipeline_dir.sh creates it and files away
# anything older runs left at the top level.
PIPELINE_DIR="/project/ANGEL/shared/Pipeline"

# Nextflow keeps its resume cache under <launch dir>/.nextflow/. Launching from
# a fixed directory is what makes -resume work across submissions: a per-job
# scratch dir would start every run from an empty cache.
RUN_DIR="${PIPELINE_DIR}/run"
mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"

nextflow run "${PIPELINE_DIR}/main.nf" \
    -resume \
    -w "${PIPELINE_DIR}/work" \
    "$@"
STATUS=$?

# Slurm profiling dumps (sps-<jobid>.tar.gz) land in the submit directory;
# file them with the job logs.
mkdir -p "${PIPELINE_DIR}/slurm/sps"
mv -f "${SLURM_SUBMIT_DIR:-$PIPELINE_DIR}"/sps-"${SLURM_JOB_ID}".tar.gz "${PIPELINE_DIR}/slurm/sps/" 2>/dev/null || true

exit $STATUS