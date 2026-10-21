// modules/create_seurat_object/main.nf

// SIF name carries the image tag, so changing params.seurat_container pulls a
// new file instead of silently reusing the old one.
def seurat_sif = "seurat_" + params.seurat_container.tokenize(':').last().replaceAll(/[^A-Za-z0-9._-]/, '_') + ".sif"

process BUILD_SEURAT_SIF {

    label 'build'
    storeDir "${projectDir}/containers"

    output:
    path("${seurat_sif}"), emit: sif
    path("${seurat_sif}.manifest"), emit: manifest

    script:
    """
    # Per-task cache/tmp: the four SIF builds run concurrently, and a shared
    # cache that each one wipes on completion deletes blobs the others are
    # still downloading. The task work dir is unique per task and per retry.
    export APPTAINER_CACHEDIR="\$PWD/.apptainer_cache"
    export APPTAINER_TMPDIR="\$PWD/.apptainer_tmp"
    mkdir -p \${APPTAINER_CACHEDIR} \${APPTAINER_TMPDIR}

    apptainer pull ${seurat_sif} ${params.seurat_container}

    # Record image digest + SIF checksum so "latest" is pinned to a known build
    bash ${projectDir}/modules/containers/record_digest.sh ${seurat_sif} "${params.seurat_container}" "\${APPTAINER_CACHEDIR}"

    # Clean up cache after successful build
    rm -rf \${APPTAINER_CACHEDIR} \${APPTAINER_TMPDIR}
    """
}


process CREATE_SEURAT_OBJECT {
    label 'build'
    tag "${library_id}"

    publishDir "${seurat_dir}/${library_id}", mode: 'copy', overwrite: true, pattern: "*.rds"

    input:
    tuple val(library_id), val(cellranger_outs), val(cellbender_dir), val(souporcell_dir), val(meta_json), val(seurat_dir)
    path(sif)

    output:
    tuple val(library_id), path("${library_id}_seurat.rds"), emit: rds
    path("${library_id}_seurat.log"),                        emit: log

    shell:
    '''
    set -euo pipefail

    # Sample-sheet metadata for this library
    cat > !{library_id}_metadata.json <<'EOF'
!{meta_json}
EOF

    {
    echo "============================================"
    echo "  CREATE_SEURAT: !{library_id}"
    echo "  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================"
    echo ""
    echo "[INFO] Container: !{sif}"
    echo ""

    apptainer exec --bind !{params.container_binds} !{sif} \
        Rscript !{projectDir}/modules/create_seurat_object/create_seurat.R \
            "!{library_id}" \
            "!{cellranger_outs}" \
            "!{cellbender_dir}" \
            "!{souporcell_dir}" \
            "!{library_id}_metadata.json"

    echo ""
    echo "[INFO] Published to: !{seurat_dir}/!{library_id}/!{library_id}_seurat.rds"
    echo "[INFO] Size: $(du -sh !{library_id}_seurat.rds | cut -f1)"
    echo "============================================"
    } 2>&1 | tee !{library_id}_seurat.log
    '''
}


process MERGE_CREATE_SEURAT_LOGS {
    label 'basic'
    tag 'merge_create_seurat_logs'

    input:
    path(log_files)
    val(logs_dir)

    output:
    path("06_Seurat.log")

    publishDir "${logs_dir}", mode: 'copy', overwrite: true

    script:
    """
    echo "================================================================" > 06_Seurat.log
    echo "  SEURAT OBJECT REPORT" >> 06_Seurat.log
    echo "  Generated: \$(date '+%Y-%m-%d %H:%M:%S')" >> 06_Seurat.log
    echo "  Samples processed: \$(ls -1 ${log_files} | wc -l)" >> 06_Seurat.log
    echo "================================================================" >> 06_Seurat.log
    echo "" >> 06_Seurat.log

    for LOG in \$(ls -1 ${log_files} | sort); do
        cat \$LOG >> 06_Seurat.log
        echo "" >> 06_Seurat.log
        echo "" >> 06_Seurat.log
    done
    """
}
