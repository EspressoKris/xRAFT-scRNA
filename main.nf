#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import groovy.json.JsonOutput

// ─── Subworkflows ───────────────────────────────────────────────────────────
include { EXTRACTION }             from './subworkflows/extraction.nf'
include { CELLRANGER_LIBRARIES }   from './subworkflows/cellranger_libraries.nf'
include { CELLRANGER }             from './subworkflows/cellranger_count.nf'
include { CELLBENDER }             from './subworkflows/cellbender.nf'
include { SOUPORCELL }             from './subworkflows/souporcell.nf'
include { CREATE_SEURAT }          from './subworkflows/create_seurat_object.nf'
include { CREATE_ANNDATA }         from './subworkflows/create_anndata_object.nf'
include { VALIDATE_OUTPUTS }       from './subworkflows/validate_outputs.nf'

// ─── Parameters ─────────────────────────────────────────────────────────────
params.samplesheet   = "${projectDir}/data/sample_sheet.csv"
params.output_format = "both"

def valid_formats = ['seurat', 'scanpy', 'both']
if (!(params.output_format in valid_formats)) {
    error "params.output_format must be one of ${valid_formats} (got '${params.output_format}')"
}

log.info """
============================================
  xRAFT-scRNA
  Started : ${new java.util.Date()}
  Output  : ${params.output_format}
  Logs    : ${params.logs_dir}
============================================
""".stripIndent()

// ─── Functions ──────────────────────────────────────────────────────────────
def parse_samplesheet(samplesheet) {
    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, strip: true)
        .filter { row -> row.library_id?.trim() }
}

// Columns that describe where data lives, not the sample itself.
// Everything else in the sample sheet is copied into the object metadata.
def path_columns = ['Raw_Data', 'fastq_dir']

// ─── Workflow ───────────────────────────────────────────────────────────────
workflow {

    ch_all = parse_samplesheet(params.samplesheet)

    // ── Extraction: only samples with Raw_Data and fastq_dir ────────────
    ch_to_extract = ch_all
        .filter { row -> row.Raw_Data?.trim() && row.fastq_dir?.trim() }
        .map { row ->
            tuple(row.library_id, row.Raw_Data.trim(), row.fastq_dir.trim())
        }

    EXTRACTION(
        ch_to_extract,
        params.logs_dir
    )

    // ── Samples that already have fastqs (no extraction needed) ─────────
    ch_existing = ch_all
        .filter { row -> !row.Raw_Data?.trim() && row.fastq_dir?.trim() }
        .map { row ->
            tuple(row.library_id, row.fastq_dir.trim())
        }

    // ── Merge: combine extracted with existing ─────────────────────────
    ch_fastqs_ready = EXTRACTION.out.fastqs.mix(ch_existing)

    ch_fastqs_ready.view { library_id, fastq_dir ->
        "Ready: ${library_id} | ${fastq_dir}"
    }

    // ── Create CellRanger Libraries CSVs ────────────────────────────────
    CELLRANGER_LIBRARIES(
        ch_fastqs_ready,
        params.cellranger_libraries,
        params.logs_dir
    )

    // ── CellRanger Count ────────────────────────────────────────────────
    ch_sampleinfo = ch_all
        .map { row ->
            tuple(row.library_id, params.transcriptome)
        }

    CELLRANGER(
        CELLRANGER_LIBRARIES.out.libraries_csv,
        ch_sampleinfo,
        params.cellranger_results,
        params.logs_dir
    )

    // ── Souporcell: genotype-based demultiplexing (skip k=1) ──────
    ch_souporcell_info = ch_all
        .filter { row -> row.souporcell_k?.trim() && row.souporcell_k.trim() as int > 1 }
        .map { row ->
            tuple(row.library_id, row.souporcell_k.trim() as int)
        }

    SOUPORCELL(
        CELLRANGER.out.outs,
        ch_souporcell_info,
        params.cellranger_results,
        params.ref_fasta,
        params.common_variants,
        params.souporcell_results,
        params.logs_dir
    )

    // ── CellBender: ambient RNA removal ───────────────────────────────
    CELLBENDER(
        CELLRANGER.out.outs,
        params.cellranger_results,
        params.cellbender_results,
        params.logs_dir
    )

    // ── Sample-sheet metadata (JSON per library) for the objects ──────
    ch_metadata = ch_all
        .map { row ->
            def meta = row.findAll { k, v -> !(k in path_columns) }
                          .collectEntries { k, v -> [(k): v?.trim() ?: ""] }
            tuple(row.library_id, JsonOutput.toJson(meta))
        }

    // ── Create Seurat Object (.rds) ───────────────────────────────────
    ch_rds        = Channel.empty()
    ch_h5ad       = Channel.empty()
    ch_seurat_sif = Channel.empty()
    ch_scanpy_sif = Channel.empty()

    if (params.output_format in ['seurat', 'both']) {
        CREATE_SEURAT(
            CELLRANGER.out.outs,
            CELLBENDER.out.filtered_h5,
            SOUPORCELL.out.results,
            ch_metadata,
            params.cellranger_results,
            params.cellbender_results,
            params.souporcell_results,
            params.seurat_results,
            params.logs_dir
        )
        ch_rds        = CREATE_SEURAT.out.rds
        ch_seurat_sif = CREATE_SEURAT.out.sif
    }

    // ── Create AnnData Object (.h5ad) ─────────────────────────────────
    if (params.output_format in ['scanpy', 'both']) {
        CREATE_ANNDATA(
            CELLRANGER.out.outs,
            CELLBENDER.out.filtered_h5,
            SOUPORCELL.out.results,
            ch_metadata,
            params.cellranger_results,
            params.cellbender_results,
            params.souporcell_results,
            params.anndata_results,
            params.logs_dir
        )
        ch_h5ad       = CREATE_ANNDATA.out.h5ad
        ch_scanpy_sif = CREATE_ANNDATA.out.sif
    }

    // ── Validate objects (structure, counts vs CellRanger, R/Python agreement)
    VALIDATE_OUTPUTS(
        ch_rds,
        ch_h5ad,
        ch_seurat_sif,
        ch_scanpy_sif,
        SOUPORCELL.out.results,
        ch_metadata,
        params.cellranger_results,
        params.souporcell_results,
        params.logs_dir
    )
}
