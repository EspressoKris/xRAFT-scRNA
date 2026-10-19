process CELLRANGER_COUNT {
    tag "${library_id}"
    label 'long'

    module 'cellranger/9.0.0'

    input:
    tuple val(library_id), path(libraries_csv), val(transcriptome), val(results_dir)

    output:
    tuple val(library_id), path("${library_id}/outs"), emit: outs
    path("${library_id}_count.log"), emit: log

    script:
    """
    set -o pipefail
    mkdir -p ${results_dir}

    {
    echo "============================================"
    echo "  CELLRANGER_COUNT: ${library_id}"
    echo "  \$(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================"
    echo ""

    if [ -d "${results_dir}/${library_id}" ]; then
        echo "[WARNING] Existing CellRanger output at ${results_dir}/${library_id} will be overwritten"
        echo ""
    fi

    echo "[INFO] Libraries CSV contents:"
    echo "--------------------------------------------"
    cat ${libraries_csv}
    echo "--------------------------------------------"
    SAMPLE_COUNT=\$(tail -n +2 ${libraries_csv} | wc -l)
    echo "[INFO] \${SAMPLE_COUNT} sample prefix(es) to process"
    echo ""

    echo "[INFO] Run parameters:"
    echo "  --id             : ${library_id}"
    echo "  --libraries      : ${libraries_csv}"
    echo "  --transcriptome  : ${transcriptome}"
    echo "  --expect-cells   : 20000"
    echo "  --create-bam     : true"
    echo "  --localcores     : ${task.cpus}"
    echo "  --localmem       : \$((${task.memory.toGiga()} - 4)) (${task.memory.toGiga()}GB allocated, 4GB headroom)"
    echo ""

    if [ ! -d "${transcriptome}" ]; then
        echo "[ERROR] Transcriptome directory not found: ${transcriptome}"
        exit 1
    fi
    echo "[INFO] Transcriptome validated: ${transcriptome}"
    echo ""

    echo "[INFO] Starting CellRanger count at \$(date)"
    echo "============================================"
    echo ""

    cellranger count \\
        --id=${library_id} \\
        --libraries=${libraries_csv} \\
        --transcriptome=${transcriptome} \\
        --expect-cells=20000 \\
        --create-bam=true \\
        --localcores=${task.cpus} \\
        --localmem=\$((${task.memory.toGiga()} - 4))

    echo ""
    echo "============================================"
    echo "  CELLRANGER_COUNT COMPLETE: ${library_id}"
    echo "============================================"
    echo "[INFO] Finished at \$(date)"
    echo ""

    if [ -f "${library_id}/outs/metrics_summary.csv" ]; then
        echo "[INFO] Metrics summary:"
        head -2 ${library_id}/outs/metrics_summary.csv | column -t -s','
        echo ""
    fi

    echo "[INFO] Output contents:"
    ls -lh ${library_id}/outs/
    echo ""

    # ── Publish results ──────────────────────────────────────────────────
    echo "[INFO] Publishing results to ${results_dir}/${library_id}..."
    rm -rf ${results_dir}/${library_id}
    mkdir -p ${results_dir}/${library_id}
    cp -r ${library_id}/outs ${results_dir}/${library_id}/
    echo "[INFO] Results published to: ${results_dir}/${library_id}"
    echo "[INFO] Size: \$(du -sh ${results_dir}/${library_id} 2>/dev/null | cut -f1)"
    echo "============================================"

    } 2>&1 | tee ${library_id}_count.log
    """
}

process MERGE_COUNT_LOGS {
    label 'basic'

    input:
    path(logs)
    val(logs_dir)

    output:
    path("03_CellRanger_count.log")

    publishDir "${logs_dir}", mode: 'copy', overwrite: true

    script:
    """
    echo "================================================================" > 03_CellRanger_count.log
    echo "  CELLRANGER COUNT REPORT" >> 03_CellRanger_count.log
    echo "  Generated: \$(date '+%Y-%m-%d %H:%M:%S')" >> 03_CellRanger_count.log
    echo "  Samples processed: \$(ls -1 *_count.log 2>/dev/null | wc -l)" >> 03_CellRanger_count.log
    echo "================================================================" >> 03_CellRanger_count.log
    echo "" >> 03_CellRanger_count.log

    for LOG in \$(ls -1 *_count.log 2>/dev/null | sort); do
        cat \$LOG >> 03_CellRanger_count.log
        echo "" >> 03_CellRanger_count.log
        echo "" >> 03_CellRanger_count.log
    done
    """
}
