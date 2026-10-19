// subworkflows/validate_outputs.nf

include { VALIDATE_SEURAT; VALIDATE_ANNDATA; COMPARE_OBJECTS; MERGE_VALIDATION_LOGS } from "${projectDir}/modules/validate/main.nf"

workflow VALIDATE_OUTPUTS {
    take:
    ch_rds          // tuple: [library_id, seurat_rds]       (empty if not built)
    ch_h5ad         // tuple: [library_id, anndata_h5ad]     (empty if not built)
    ch_seurat_sif   // path: seurat SIF                      (empty if not built)
    ch_scanpy_sif   // path: scanpy SIF                      (empty if not built)
    ch_souporcell   // tuple: [library_id, souporcell_output_dir]
    ch_metadata     // tuple: [library_id, metadata_json]
    cellranger_dir  // val: published CellRanger results
    souporcell_dir  // val: published Souporcell results
    logs_dir        // val: logs directory

    main:
    // Published Souporcell dir, or "NONE" for k=1 libraries
    ch_sc = ch_souporcell.map { library_id, sc -> tuple(library_id, "${souporcell_dir}/${library_id}") }

    // [library_id, object, cellranger_outs, souporcell_dir|NONE, meta_json]
    def with_context = { ch_obj ->
        ch_obj
            .join(ch_sc, remainder: true)
            .filter { it[1] != null }                       // drop Souporcell-only rows
            .map { library_id, obj, sc -> tuple(library_id, obj, sc ?: "NONE") }
            .join(ch_metadata)
            .map { library_id, obj, sc, meta ->
                tuple(library_id, obj, "${cellranger_dir}/${library_id}/outs", sc, meta)
            }
    }

    VALIDATE_SEURAT(with_context(ch_rds), ch_seurat_sif)
    VALIDATE_ANNDATA(with_context(ch_h5ad), ch_scanpy_sif)

    // Only fires for libraries that have both objects
    COMPARE_OBJECTS(VALIDATE_SEURAT.out.summary.join(VALIDATE_ANNDATA.out.summary))

    MERGE_VALIDATION_LOGS(
        VALIDATE_SEURAT.out.log
            .mix(VALIDATE_ANNDATA.out.log)
            .mix(COMPARE_OBJECTS.out.log)
            .collect(),
        logs_dir
    )

    emit:
    seurat_summary  = VALIDATE_SEURAT.out.summary
    anndata_summary = VALIDATE_ANNDATA.out.summary
}
