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

# Nextflow keeps its resume cache (an embedded LevelDB) under <launch dir>/
# .nextflow/. Launching from a FIXED directory is what makes -resume work
# across submissions. The DB needs a file lock; if the filesystem refuses one
# ("Can't open cache DB"), the broken session cache is moved aside and the run
# starts with a fresh cache (finished work is recomputed once).
RUN_DIR="${PIPELINE_DIR}/run"
mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"

# A previous coordinator killed by scancel can leave a stale LevelDB lock.
# Only this job runs Nextflow against this cache, so clear it.
if [ "$(squeue -h -u "$USER" -n xraft_scrna -o %i | grep -vc "^${SLURM_JOB_ID}$")" -eq 0 ]; then
    rm -f "${RUN_DIR}"/.nextflow/cache/*/db/LOCK
fi

# ── Run Nextflow under a watchdog ────────────────────────────────────────────
# If the Slurm task monitor stops polling (no "Task monitor" log line for
# STALL_MINUTES while the process is still alive), kill Nextflow and relaunch
# with -resume. Finished tasks are cached, so nothing is recomputed.
STALL_MINUTES=90
MAX_RELAUNCH=5
NF_LOG="${RUN_DIR}/.nextflow.log"
STATUS=1

for attempt in $(seq 0 ${MAX_RELAUNCH}); do
    [ "$attempt" -gt 0 ] && echo "[WATCHDOG] relaunch ${attempt}/${MAX_RELAUNCH} at $(date '+%F %T')"
    nextflow run "${PIPELINE_DIR}/main.nf" -resume -w "${PIPELINE_DIR}/work" "$@" &
    NF_PID=$!

    while kill -0 "$NF_PID" 2>/dev/null; do
        sleep 300
        [ -f "$NF_LOG" ] || continue
        last_poll=$(grep -E "\[Task monitor\]" "$NF_LOG" | tail -1 | cut -c1-19)
        [ -n "$last_poll" ] || continue
        last_epoch=$(date -d "${last_poll/-/ } $(date +%Y)" +%s 2>/dev/null) || continue
        age=$(( ($(date +%s) - last_epoch) / 60 ))
        if [ "$age" -ge "$STALL_MINUTES" ]; then
            echo "[WATCHDOG] task monitor silent for ${age} min (last: ${last_poll}); killing Nextflow pid ${NF_PID}"
            kill "$NF_PID" 2>/dev/null; sleep 30; kill -9 "$NF_PID" 2>/dev/null
            wait "$NF_PID" 2>/dev/null
            # LevelDB lock on the cache DB can outlive the process on Ceph;
            # the PID is confirmed dead, so clear it before relaunching.
            sleep 10; rm -f "${RUN_DIR}"/.nextflow/cache/*/db/LOCK
            STATUS=75   # EX_TEMPFAIL -> relaunch
            break
        fi
    done

    if [ "$STATUS" -ne 75 ]; then
        wait "$NF_PID"; STATUS=$?
        # Cache DB could not be locked/opened: park that session's cache and
        # relaunch once with a fresh one rather than failing the whole run.
        if [ "$STATUS" -ne 0 ] && grep -q "Can't open cache DB" "$NF_LOG" 2>/dev/null \
           && [ "${CACHE_RESET:-0}" -eq 0 ]; then
            CACHE_RESET=1
            echo "[WATCHDOG] cache DB unusable; moving .nextflow/cache aside and relaunching with a fresh cache (finished tasks will be recomputed)"
            mv "${RUN_DIR}/.nextflow/cache" "${RUN_DIR}/.nextflow/cache.broken.$(date +%Y%m%d_%H%M%S)" 2>/dev/null
            rm -f "${RUN_DIR}/.nextflow/history"
            continue
        fi
        break
    fi
    STATUS=1
done

# Slurm profiling dumps (sps-<jobid>.tar.gz) land in the submit directory;
# file them with the job logs.
mkdir -p "${PIPELINE_DIR}/slurm/sps"
mv -f "${SLURM_SUBMIT_DIR:-$PIPELINE_DIR}"/sps-"${SLURM_JOB_ID}".tar.gz "${PIPELINE_DIR}/slurm/sps/" 2>/dev/null || true

exit $STATUS
