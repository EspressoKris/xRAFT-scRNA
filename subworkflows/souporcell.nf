include { BUILD_SOUPORCELL_SIF; SOUPORCELL_DEMUX; MERGE_SOUPORCELL_LOGS } from "${projectDir}/modules/souporcell/main.nf"

workflow SOUPORCELL {
    take:
    ch_cellranger_outs  // tuple: [library_id, cellranger_outs_dir]
    ch_souporcell_info  // tuple: [library_id, souporcell_k]
    results_dir         // val: path to CellRanger results (for BAM/barcodes)
    ref_fasta           // val: path to reference genome fasta
    common_variants     // val: path to common-variants VCF
    souporcell_dir      // val: path to Souporcell output directory
    logs_dir            // val: path to logs directory

    main:
    // Build (or reuse cached) Souporcell SIF
    BUILD_SOUPORCELL_SIF()

    // Wrap val inputs as value channels for safe use in combine
    ch_results_dir     = Channel.value(results_dir)
    ch_souporcell_dir  = Channel.value(souporcell_dir)
    ch_ref_fasta       = Channel.value(ref_fasta)
    ch_common_variants = Channel.value(common_variants)

    // Use the PUBLISHED results path for BAM/barcodes
    // Join with souporcell_k from sample sheet
    ch_input = ch_cellranger_outs
        .combine(ch_results_dir)
        .map { library_id, outs_dir, res_dir ->
            tuple(library_id, "${res_dir}/${library_id}/outs")
        }
        .join(ch_souporcell_info)
        .combine(ch_souporcell_dir)
        .combine(ch_ref_fasta)
        .combine(ch_common_variants)
        .map { library_id, outs_dir, souporcell_k, sc_dir, ref, vcf ->
            tuple(library_id, outs_dir, sc_dir, souporcell_k, ref, vcf)
        }

    SOUPORCELL_DEMUX(ch_input, BUILD_SOUPORCELL_SIF.out.sif)

    MERGE_SOUPORCELL_LOGS(
        SOUPORCELL_DEMUX.out.log.collect(),
        logs_dir
    )

    emit:
    results = SOUPORCELL_DEMUX.out.results  // tuple: [library_id, souporcell_output_dir]
}
