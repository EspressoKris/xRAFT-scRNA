// subworkflows/create_seurat_object.nf

include { BUILD_SEURAT_SIF; CREATE_SEURAT_OBJECT; MERGE_CREATE_SEURAT_LOGS } from "${projectDir}/modules/create_seurat_object/main.nf"

workflow CREATE_SEURAT {
    take:
    ch_cellranger_outs  // tuple: [library_id, cellranger_outs_dir]
    ch_cellbender_h5    // tuple: [library_id, cellbender_filtered_h5]
    ch_souporcell       // tuple: [library_id, souporcell_output_dir]
    ch_metadata         // tuple: [library_id, metadata_json]  (sample-sheet columns)
    cellranger_dir      // val: path to published CellRanger results
    cellbender_dir      // val: path to published CellBender results
    souporcell_dir      // val: path to published Souporcell results
    seurat_dir          // val: path to Seurat output directory
    logs_dir            // val: path to logs directory

    main:
    // Build (or reuse cached) Seurat SIF
    BUILD_SEURAT_SIF()

    // Use PUBLISHED paths rather than Nextflow work dirs
    ch_cr = ch_cellranger_outs
        .map { library_id, outs_dir ->
            tuple(library_id, "${cellranger_dir}/${library_id}/outs")
        }

    ch_cb = ch_cellbender_h5
        .map { library_id, h5_file ->
            tuple(library_id, "${cellbender_dir}/${library_id}")
        }

    ch_sc = ch_souporcell
        .map { library_id, sc_dir ->
            tuple(library_id, "${souporcell_dir}/${library_id}")
        }

    // Join CellRanger + CellBender (both always present)
    // Then left-join Souporcell: remainder = true so samples without
    // Souporcell (k=1) still proceed with souporcell_dir = "NONE"
    ch_input = ch_cr
        .join(ch_cb)
        .join(ch_sc, remainder: true)
        .map { library_id, cr_outs, cb_dir, sc_dir ->
            tuple(library_id, cr_outs, cb_dir, sc_dir ?: "NONE")
        }
        .join(ch_metadata)
        .map { library_id, cr_outs, cb_dir, sc_dir, meta_json ->
            tuple(library_id, cr_outs, cb_dir, sc_dir, meta_json, seurat_dir)
        }

    CREATE_SEURAT_OBJECT(ch_input, BUILD_SEURAT_SIF.out.sif)

    MERGE_CREATE_SEURAT_LOGS(
        CREATE_SEURAT_OBJECT.out.log.collect(),
        logs_dir
    )

    emit:
    rds = CREATE_SEURAT_OBJECT.out.rds  // tuple: [library_id, seurat_rds]
    sif  = BUILD_SEURAT_SIF.out.sif
}
