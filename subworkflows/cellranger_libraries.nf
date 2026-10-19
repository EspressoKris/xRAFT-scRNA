include { CREATE_CELLRANGER_LIBRARIES; MERGE_LIBRARIES_LOGS } from "${projectDir}/modules/cellranger_libraries/main.nf"

workflow CELLRANGER_LIBRARIES {
    take:
    ch_fastqs      // tuple: [library_id, fastq_dir]
    libraries_dir  // val: path to CellRanger_Libraries folder
    logs_dir       // val: path to logs directory

    main:
    ch_input = ch_fastqs.map { library_id, fastq_dir ->
        tuple(library_id, fastq_dir, libraries_dir)
    }

    CREATE_CELLRANGER_LIBRARIES(ch_input)

    MERGE_LIBRARIES_LOGS(
        CREATE_CELLRANGER_LIBRARIES.out.log.collect(),
        logs_dir
    )

    emit:
    libraries_csv = CREATE_CELLRANGER_LIBRARIES.out.libraries_csv
}
