// SIF name carries the image tag, so changing params.souporcell_container pulls a
// new file instead of silently reusing the old one.
def souporcell_sif = "souporcell_" + params.souporcell_container.tokenize(':').last().replaceAll(/[^A-Za-z0-9._-]/, '_') + ".sif"

process BUILD_SOUPORCELL_SIF {

    label 'build'
    storeDir "${projectDir}/containers"

    output:
    path("${souporcell_sif}"), emit: sif
    path("${souporcell_sif}.manifest"), emit: manifest

    script:
    """
    # Per-task cache/tmp: the four SIF builds run concurrently, and a shared
    # cache that each one wipes on completion deletes blobs the others are
    # still downloading. The task work dir is unique per task and per retry.
    export APPTAINER_CACHEDIR="\$PWD/.apptainer_cache"
    export APPTAINER_TMPDIR="\$PWD/.apptainer_tmp"
    mkdir -p \${APPTAINER_CACHEDIR} \${APPTAINER_TMPDIR}

    singularity pull ${souporcell_sif} ${params.souporcell_container}

    # Record image digest + SIF checksum so "latest" is pinned to a known build
    bash ${projectDir}/modules/containers/record_digest.sh ${souporcell_sif} "${params.souporcell_container}" "\${APPTAINER_CACHEDIR}"

    # Clean up cache after successful build
    rm -rf \${APPTAINER_CACHEDIR} \${APPTAINER_TMPDIR}
    """
}

process SOUPORCELL_DEMUX {

    tag "${library_id}"
    label 'long'

    input:
    tuple val(library_id), val(outs_dir), val(souporcell_dir), val(souporcell_k), val(ref_fasta), val(common_variants)
    path(sif)

    output:
    tuple val(library_id), path("${library_id}_souporcell/"), emit: results
    path("${library_id}_souporcell.log"),                     emit: log

    script:
    def output_dir  = "${souporcell_dir}/${library_id}"
    def bam         = "${outs_dir}/possorted_genome_bam.bam"
    def barcodes    = "${outs_dir}/filtered_feature_bc_matrix/barcodes.tsv.gz"
    """
    set -euo pipefail

    # ── Validate inputs ──────────────────────────────────────────────────
    if [ ! -f "${bam}" ]; then
        echo "[ERROR] BAM file not found: ${bam}" | tee ${library_id}_souporcell.log
        echo "[ERROR] CellRanger outs dir: ${outs_dir}" | tee -a ${library_id}_souporcell.log
        ls -la ${outs_dir}/ 2>&1 | tee -a ${library_id}_souporcell.log
        exit 1
    fi

    if [ ! -f "${barcodes}" ]; then
        echo "[ERROR] Barcodes file not found: ${barcodes}" | tee ${library_id}_souporcell.log
        exit 1
    fi

    if [ ! -f "${common_variants}" ]; then
        echo "[ERROR] Common-variants VCF not found: ${common_variants}" | tee ${library_id}_souporcell.log
        exit 1
    fi

    # ── Full run ─────────────────────────────────────────────────────────
    {
    echo "============================================"
    echo "  SOUPORCELL: ${library_id}"
    echo "  \$(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================"
    echo ""

    if [ -d "${output_dir}" ]; then
        echo "[WARNING] Existing Souporcell output at ${output_dir} will be overwritten"
        echo ""
    fi

    echo "[INFO] Run parameters:"
    echo "  -i (BAM)         : ${bam}"
    echo "  -b (barcodes)    : ${barcodes}"
    echo "  -f (ref fasta)   : ${ref_fasta}"
    echo "  -k (clusters)    : ${souporcell_k}"
    echo "  -t (threads)     : ${task.cpus}"
    echo "  -o (output)      : ${library_id}_souporcell"
    echo "  --common_variants: ${common_variants}"
    echo "  --skip_remap     : SKIP_REMAP"
    echo "  Container        : ${sif}"
    echo ""

    echo "[INFO] Starting Souporcell at \$(date)"
    echo "============================================"
    echo ""

    singularity exec --bind ${params.container_binds} ${sif} souporcell_pipeline.py \\
        -i ${bam} \\
        -b ${barcodes} \\
        -f ${ref_fasta} \\
        -t ${task.cpus} \\
        -o ${library_id}_souporcell \\
        --common_variants ${common_variants} \\
        -k ${souporcell_k} \\
        --skip_remap SKIP_REMAP

    echo ""
    echo "============================================"
    echo "  SOUPORCELL COMPLETE: ${library_id}"
    echo "============================================"
    echo "[INFO] Finished at \$(date)"
    echo ""

    echo "[INFO] Output files:"
    ls -lh ${library_id}_souporcell/
    echo ""

    # Delete intermediate .bam and .bai to save space
    find ${library_id}_souporcell -type f \\( -name "*.bam" -o -name "*.bai" \\) -exec rm -f {} \\;
    echo "[INFO] Cleaned up intermediate BAM/BAI files"
    echo ""

    # Publish results
    rm -rf ${output_dir}
    mkdir -p ${output_dir}
    cp -r ${library_id}_souporcell/* ${output_dir}/

    echo "[INFO] Results published to: ${output_dir}"
    echo "[INFO] Total size: \$(du -sh ${output_dir} 2>/dev/null | cut -f1)"
    echo "============================================"

    } 2>&1 | tee ${library_id}_souporcell.log
    """
}

process MERGE_SOUPORCELL_LOGS {

    label 'basic'

    input:
    path(log_files)
    val(logs_dir)

    output:
    path("05_Souporcell.log")

    publishDir "${logs_dir}", mode: 'copy', overwrite: true

    script:
    """
    echo "================================================================" > 05_Souporcell.log
    echo "  SOUPORCELL REPORT" >> 05_Souporcell.log
    echo "  Generated: \$(date '+%Y-%m-%d %H:%M:%S')" >> 05_Souporcell.log
    echo "  Samples processed: \$(ls -1 *_souporcell.log 2>/dev/null | wc -l)" >> 05_Souporcell.log
    echo "================================================================" >> 05_Souporcell.log
    echo "" >> 05_Souporcell.log

    for LOG in \$(ls -1 *_souporcell.log 2>/dev/null | sort); do
        cat \$LOG >> 05_Souporcell.log
        echo "" >> 05_Souporcell.log
        echo "" >> 05_Souporcell.log
    done
    """
}
