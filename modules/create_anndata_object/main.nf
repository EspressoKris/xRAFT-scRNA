// modules/create_anndata_object/main.nf

// SIF name carries the image tag, so changing params.scanpy_container pulls a
// new file instead of silently reusing the old one.
def scanpy_sif = "scanpy_" + params.scanpy_container.tokenize(':').last().replaceAll(/[^A-Za-z0-9._-]/, '_') + ".sif"

process BUILD_SCANPY_SIF {

    label 'build'
    storeDir "${projectDir}/containers"

    output:
    path("${scanpy_sif}"), emit: sif
    path("${scanpy_sif}.manifest"), emit: manifest

    script:
    """
    export APPTAINER_CACHEDIR="${projectDir}/containers/.apptainer_cache"
    export APPTAINER_TMPDIR="${projectDir}/containers/.apptainer_tmp"
    mkdir -p \${APPTAINER_CACHEDIR} \${APPTAINER_TMPDIR}

    apptainer pull ${scanpy_sif} ${params.scanpy_container}

    # Record image digest + SIF checksum so "latest" is pinned to a known build
    bash ${projectDir}/modules/containers/record_digest.sh ${scanpy_sif} "${params.scanpy_container}" "\${APPTAINER_CACHEDIR}"

    # Clean up cache after successful build
    rm -rf \${APPTAINER_CACHEDIR} \${APPTAINER_TMPDIR}
    """
}


process CREATE_ANNDATA_OBJECT {
    label 'build'
    tag "${library_id}"

    publishDir "${anndata_dir}/${library_id}", mode: 'copy', overwrite: true, pattern: "*.h5ad"

    input:
    tuple val(library_id), val(cellranger_outs), val(cellbender_dir), val(souporcell_dir), val(meta_json), val(anndata_dir)
    path(sif)

    output:
    tuple val(library_id), path("${library_id}.h5ad"), emit: h5ad
    path("${library_id}_anndata.log"),                 emit: log

    shell:
    '''
    set -euo pipefail

    # Sample-sheet metadata for this library
    cat > !{library_id}_metadata.json <<'EOF'
!{meta_json}
EOF

    {
    echo "============================================"
    echo "  CREATE_ANNDATA: !{library_id}"
    echo "  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================"
    echo ""
    echo "[INFO] Container: !{sif}"
    echo ""

    apptainer exec --bind !{params.container_binds} !{sif} \
        python3 !{projectDir}/modules/create_anndata_object/create_anndata.py \
            "!{library_id}" \
            "!{cellranger_outs}" \
            "!{cellbender_dir}" \
            "!{souporcell_dir}" \
            "!{library_id}_metadata.json"

    echo ""
    echo "[INFO] Published to: !{anndata_dir}/!{library_id}/!{library_id}.h5ad"
    echo "[INFO] Size: $(du -sh !{library_id}.h5ad | cut -f1)"
    echo "============================================"
    } 2>&1 | tee !{library_id}_anndata.log
    '''
}


process MERGE_CREATE_ANNDATA_LOGS {
    label 'basic'
    tag 'merge_create_anndata_logs'

    input:
    path(log_files)
    val(logs_dir)

    output:
    path("07_AnnData.log")

    publishDir "${logs_dir}", mode: 'copy', overwrite: true

    script:
    """
    echo "================================================================" > 07_AnnData.log
    echo "  ANNDATA OBJECT REPORT" >> 07_AnnData.log
    echo "  Generated: \$(date '+%Y-%m-%d %H:%M:%S')" >> 07_AnnData.log
    echo "  Samples processed: \$(ls -1 *_anndata.log 2>/dev/null | wc -l)" >> 07_AnnData.log
    echo "================================================================" >> 07_AnnData.log
    echo "" >> 07_AnnData.log

    for LOG in \$(ls -1 *_anndata.log 2>/dev/null | sort); do
        cat \$LOG >> 07_AnnData.log
        echo "" >> 07_AnnData.log
        echo "" >> 07_AnnData.log
    done
    """
}
