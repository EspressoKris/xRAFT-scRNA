// modules/validate/main.nf

process VALIDATE_SEURAT {
    label 'build'
    tag "${library_id}"

    input:
    tuple val(library_id), path(rds), val(cellranger_outs), val(souporcell_dir), val(meta_json)
    path(sif)

    output:
    tuple val(library_id), path("${library_id}_seurat_summary.json"), emit: summary
    path("${library_id}_validate_seurat.log"),                        emit: log

    shell:
    '''
    set -euo pipefail

    cat > !{library_id}_metadata.json <<'EOF'
!{meta_json}
EOF

    {
    echo "============================================"
    echo "  VALIDATE_SEURAT: !{library_id}"
    echo "  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================"

    apptainer exec --bind !{params.container_binds} !{sif} \
        Rscript !{projectDir}/modules/validate/validate_seurat.R \
            "!{library_id}" "!{rds}" "!{cellranger_outs}" "!{souporcell_dir}" "!{library_id}_metadata.json"

    echo "============================================"
    } 2>&1 | tee !{library_id}_validate_seurat.log
    '''
}


process VALIDATE_ANNDATA {
    label 'build'
    tag "${library_id}"

    input:
    tuple val(library_id), path(h5ad), val(cellranger_outs), val(souporcell_dir), val(meta_json)
    path(sif)

    output:
    tuple val(library_id), path("${library_id}_anndata_summary.json"), emit: summary
    path("${library_id}_validate_anndata.log"),                        emit: log

    shell:
    '''
    set -euo pipefail

    cat > !{library_id}_metadata.json <<'EOF'
!{meta_json}
EOF

    {
    echo "============================================"
    echo "  VALIDATE_ANNDATA: !{library_id}"
    echo "  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================"

    apptainer exec --bind !{params.container_binds} !{sif} \
        python3 !{projectDir}/modules/validate/validate_anndata.py \
            "!{library_id}" "!{h5ad}" "!{cellranger_outs}" "!{souporcell_dir}" "!{library_id}_metadata.json"

    echo "============================================"
    } 2>&1 | tee !{library_id}_validate_anndata.log
    '''
}


process COMPARE_OBJECTS {
    label 'basic'
    tag "${library_id}"

    input:
    tuple val(library_id), path(seurat_summary), path(anndata_summary)

    output:
    path("${library_id}_compare.log"), emit: log

    script:
    """
    set -euo pipefail
    {
    echo "============================================"
    echo "  COMPARE_OBJECTS: ${library_id}"
    echo "  \$(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================"

    python3 ${projectDir}/modules/validate/compare_summaries.py \\
        "${library_id}" "${seurat_summary}" "${anndata_summary}"

    echo "============================================"
    } 2>&1 | tee ${library_id}_compare.log
    """
}


process MERGE_VALIDATION_LOGS {
    label 'basic'
    tag 'merge_validation_logs'

    input:
    path(log_files)
    val(logs_dir)

    output:
    path("08_Validation.log")

    publishDir "${logs_dir}", mode: 'copy', overwrite: true

    script:
    """
    echo "================================================================" > 08_Validation.log
    echo "  OBJECT VALIDATION REPORT" >> 08_Validation.log
    echo "  Generated: \$(date '+%Y-%m-%d %H:%M:%S')" >> 08_Validation.log
    echo "  Checks run: \$(ls -1 *.log 2>/dev/null | grep -v '^08_' | wc -l)" >> 08_Validation.log
    echo "  Failures  : \$(grep -l '\\[FAIL\\]' *_validate_*.log *_compare.log 2>/dev/null | wc -l)" >> 08_Validation.log
    echo "================================================================" >> 08_Validation.log
    echo "" >> 08_Validation.log

    for LOG in \$(ls -1 *_validate_seurat.log *_validate_anndata.log *_compare.log 2>/dev/null | sort); do
        cat \$LOG >> 08_Validation.log
        echo "" >> 08_Validation.log
        echo "" >> 08_Validation.log
    done
    """
}
