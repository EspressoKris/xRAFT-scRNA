// SIF name carries the image tag, so changing params.cellbender_container pulls a
// new file instead of silently reusing the old one.
def cellbender_sif = "cellbender_" + params.cellbender_container.tokenize(':').last().replaceAll(/[^A-Za-z0-9._-]/, '_') + ".sif"

process BUILD_CELLBENDER_SIF {

    label 'build'
    storeDir "${projectDir}/containers"

    output:
    path("${cellbender_sif}"), emit: sif
    path("${cellbender_sif}.manifest"), emit: manifest

    script:
    """
    # Per-task cache/tmp: the four SIF builds run concurrently, and a shared
    # cache that each one wipes on completion deletes blobs the others are
    # still downloading. The task work dir is unique per task and per retry.
    export APPTAINER_CACHEDIR="\$PWD/.apptainer_cache"
    export APPTAINER_TMPDIR="\$PWD/.apptainer_tmp"
    mkdir -p \${APPTAINER_CACHEDIR} \${APPTAINER_TMPDIR}

    apptainer pull ${cellbender_sif} ${params.cellbender_container}

    # Record image digest + SIF checksum so "latest" is pinned to a known build
    bash ${projectDir}/modules/containers/record_digest.sh ${cellbender_sif} "${params.cellbender_container}" "\${APPTAINER_CACHEDIR}"

    # Clean up cache after successful build
    rm -rf \${APPTAINER_CACHEDIR} \${APPTAINER_TMPDIR}
    """
}

process CELLBENDER_REMOVE_BACKGROUND {

    tag "${library_id}"
    label 'gpu'

    input:
    tuple val(library_id), val(outs_dir), val(results_dir), val(cellbender_dir)
    path(sif)

    output:
    tuple val(library_id), path("${library_id}_cellbender_output.h5"),          emit: h5
    tuple val(library_id), path("${library_id}_cellbender_output_filtered.h5"), emit: filtered_h5
    path("${library_id}_cellbender.log"),                                       emit: log

    script:
    def output_dir = "${cellbender_dir}/${library_id}"
    def raw_h5     = "${outs_dir}/raw_feature_bc_matrix.h5"
    """
    set -euo pipefail

    # ── Validate input ───────────────────────────────────────────────────
    if [ ! -f "${raw_h5}" ]; then
        echo "[ERROR] Input file not found: ${raw_h5}" | tee ${library_id}_cellbender.log
        echo "[ERROR] CellRanger outs dir: ${outs_dir}" | tee -a ${library_id}_cellbender.log
        ls -la ${outs_dir}/ 2>&1 | tee -a ${library_id}_cellbender.log
        exit 1
    fi

    # ── Full run ─────────────────────────────────────────────────────────
    {
    echo "============================================"
    echo "  CELLBENDER: ${library_id}"
    echo "  \$(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================"
    echo ""

    if [ -d "${output_dir}" ]; then
        echo "[WARNING] Existing CellBender output at ${output_dir} will be overwritten"
        echo ""
    fi

    echo "[INFO] Run parameters:"
    echo "  --input          : ${raw_h5}"
    echo "  --output         : ${library_id}_cellbender_output.h5"
    echo "  --cuda           : true"
    echo "  Container        : ${sif}"
    echo ""

    echo "[INFO] Starting CellBender remove-background at \$(date)"
    echo "============================================"
    echo ""

    apptainer exec --nv --writable-tmpfs --bind ${params.container_binds} ${sif} \\
        cellbender remove-background \\
            --cuda \\
            --input ${raw_h5} \\
            --output ${library_id}_cellbender_output.h5

    echo ""
    echo "============================================"
    echo "  CELLBENDER COMPLETE: ${library_id}"
    echo "============================================"
    echo "[INFO] Finished at \$(date)"
    echo ""

    echo "[INFO] Output files:"
    ls -lh ${library_id}_cellbender_output*.h5
    echo ""

    # Publish results
    rm -rf ${output_dir}
    mkdir -p ${output_dir}
    cp ${library_id}_cellbender_output*.h5 ${output_dir}/

    echo "[INFO] Results published to: ${output_dir}"
    echo "[INFO] Total size: \$(du -sh ${output_dir} 2>/dev/null | cut -f1)"
    echo "============================================"

    } 2>&1 | tee ${library_id}_cellbender.log
    """
}

process MERGE_CELLBENDER_LOGS {

    label 'basic'

    input:
    path(log_files)
    val(logs_dir)

    output:
    path("04_CellBender.log")

    publishDir "${logs_dir}", mode: 'copy', overwrite: true

    script:
    """
    echo "================================================================" > 04_CellBender.log
    echo "  CELLBENDER REPORT" >> 04_CellBender.log
    echo "  Generated: \$(date '+%Y-%m-%d %H:%M:%S')" >> 04_CellBender.log
    echo "  Samples processed: \$(ls -1 ${log_files} | wc -l)" >> 04_CellBender.log
    echo "================================================================" >> 04_CellBender.log
    echo "" >> 04_CellBender.log

    for LOG in \$(ls -1 ${log_files} | sort); do
        cat \$LOG >> 04_CellBender.log
        echo "" >> 04_CellBender.log
        echo "" >> 04_CellBender.log
    done
    """
}
