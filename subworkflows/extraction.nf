include { EXTRACT_RAW_DATA; MERGE_EXTRACT_LOGS } from "${projectDir}/modules/extract/main.nf"

workflow EXTRACTION {
    take:
    ch_samples      // tuple: [library_id, tar_files, fastq_dir]
    logs_dir        // val: path to logs directory

    main:
    EXTRACT_RAW_DATA(ch_samples)

    MERGE_EXTRACT_LOGS(
        EXTRACT_RAW_DATA.out.log.collect(),
        logs_dir
    )

    emit:
    fastqs = EXTRACT_RAW_DATA.out.fastqs  // tuple: [library_id, fastq_dir]
}
