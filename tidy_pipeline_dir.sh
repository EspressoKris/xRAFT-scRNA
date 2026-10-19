#!/bin/bash
# ─── Tidy the pipeline directory ─────────────────────────────────────────────
# Files away everything Slurm and Nextflow leave at the top level so the
# directory only holds the pipeline itself:
#
#   Report-*.out / *.out, *.err  →  slurm/
#   sps-*.tar.gz, sps-*/ (profiling) →  slurm/sps/
#   .nextflow.log*               →  run/
#   containers/*.sif without a .manifest get one (digest unknown, checksum recorded)
#
# Safe to run any time, including while a pipeline job is queued (it never
# touches work/, containers/, resources/, logs/ or reports/).

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

mkdir -p slurm/sps run

moved=0
shopt -s nullglob
for f in Report-*.out Report-*.err *.out *.err; do
    [ -f "$f" ] || continue
    mv -f "$f" slurm/ && moved=$((moved + 1))
done
for f in sps-*.tar.gz sps-*/; do
    [ -e "$f" ] || continue
    mv -f "${f%/}" slurm/sps/ && moved=$((moved + 1))
done
for f in .nextflow.log*; do
    mv -f "$f" run/ && moved=$((moved + 1))
done

# ── Manifests for SIFs pulled before manifests existed ───────────────────────
# Without one, Nextflow's storeDir would re-pull the image on the next run.
# The image URI is taken from nextflow.config when the SIF name matches the
# pinned tag (e.g. seurat_5.5.1.sif <-> seurat_container ...:5.5.1), else unknown.
for sif in containers/*.sif; do
    [ -f "$sif" ] && [ ! -f "${sif}.manifest" ] || continue
    base=$(basename "$sif" .sif); key="${base%%_*}"
    img=$(grep -E "^\s*${key}_container\s*=" nextflow.config 2>/dev/null | sed -E 's/.*"([^"]+)".*/\1/')
    tag="${img##*:}"
    if [ -n "$img" ] && [ "${base#${key}_}" = "${tag//[^A-Za-z0-9._-]/_}" ]; then
        bash modules/containers/record_digest.sh "$sif" "$img" --existing
    else
        bash modules/containers/record_digest.sh "$sif" "unknown" --existing
    fi
done

echo "Moved ${moved} file(s)."
echo "  slurm/     : $(ls -1 slurm 2>/dev/null | grep -vc '^sps$' || true) job logs"
echo "  slurm/sps/ : $(ls -1 slurm/sps 2>/dev/null | wc -l) profiling tarballs"
echo "  run/       : $(ls -1A run 2>/dev/null | grep -c nextflow.log || true) nextflow logs"
