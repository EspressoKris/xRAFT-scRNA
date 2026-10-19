#!/usr/bin/env Rscript

# ─── Validate an xRAFT-scRNA Seurat object against its inputs ───────────────
# Mirror of validate_anndata.py. Checks (all must pass, exit 1 otherwise):
#   structure  : single RNA Assay5 with layers counts / RNA_raw_counts / RNA_cellbender_counts
#   identity   : counts == RNA_raw_counts exactly
#   values     : counts finite, non-negative, integer-valued
#   cellbender : 0 <= RNA_cellbender_counts <= RNA_raw_counts element-wise
#   barcodes   : colnames == CellRanger filtered barcodes.tsv.gz (same order)
#   features   : rownames == CellRanger features (made unique, same order)
#   counts     : sum(RNA_raw_counts) == sum(CellRanger matrix.mtx.gz)
#   souporcell : status/assignment present; k>1 => matches clusters.tsv; k==1 => all NA
#   metadata   : every sample-sheet column present with exactly one value
# Writes <library_id>_seurat_summary.json for cross-format comparison.

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
    stop("Usage: validate_seurat.R <library_id> <rds> <cellranger_outs> <souporcell_dir|NONE> <metadata_json>")
}
library_id      <- args[1]
rds_path        <- args[2]
cellranger_outs <- args[3]
souporcell_dir  <- args[4]
metadata_json   <- args[5]
mtx_dir         <- file.path(cellranger_outs, "filtered_feature_bc_matrix")

VALID_STATUS <- c("singlet", "doublet", "unassigned")
failures <- character(0)

check <- function(cond, msg) {
    ok <- isTRUE(cond)
    cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", msg))
    if (!ok) failures <<- c(failures, msg)
    invisible(ok)
}
md5_of <- function(x) {
    # match python: hashlib.md5("\n".join(items)) -- no trailing newline
    f <- tempfile()
    cat(paste(as.character(x), collapse = "\n"), file = f)
    unname(tools::md5sum(f))
}

cat("=== Validating Seurat object:", library_id, "===\n")
cat("rds             :", rds_path, "\n")
cat("CellRanger outs :", cellranger_outs, "\n")
cat("Souporcell dir  :", souporcell_dir, "\n")

# ─── Load ───────────────────────────────────────────────────────────────────
obj <- readRDS(rds_path)
cat("\nLoaded:", ncol(obj), "cells x", nrow(obj), "features\n")

# ─── Structure ──────────────────────────────────────────────────────────────
cat("\n[structure]\n")
check(identical(Assays(obj), "RNA"), "single assay named RNA")
check(inherits(obj[["RNA"]], "Assay5"), "RNA assay is Assay5")
layers <- Layers(obj[["RNA"]])
for (l in c("counts", "RNA_raw_counts", "RNA_cellbender_counts")) {
    check(l %in% layers, sprintf("layer '%s' present", l))
}
if (length(failures) > 0) stop("Missing layers; aborting further checks.")

counts <- as(LayerData(obj, assay = "RNA", layer = "counts"), "CsparseMatrix")
raw    <- as(LayerData(obj, assay = "RNA", layer = "RNA_raw_counts"), "CsparseMatrix")
cb     <- as(LayerData(obj, assay = "RNA", layer = "RNA_cellbender_counts"), "CsparseMatrix")

check(identical(dim(counts), dim(raw)) && identical(dim(raw), dim(cb)), "all layers share dimensions")
check(identical(dimnames(raw), dimnames(cb)) && identical(dimnames(counts), dimnames(raw)), "all layers share dimnames")
check(ncol(obj) > 0 && nrow(obj) > 0, "non-empty object")

# ─── Identity / values ──────────────────────────────────────────────────────
cat("\n[values]\n")
check(all(counts == raw), "counts is identical to RNA_raw_counts")
check(all(is.finite(raw@x)) && all(is.finite(cb@x)), "all counts finite")
check(all(raw@x >= 0) && all(cb@x >= 0), "all counts non-negative")
check(isTRUE(all.equal(raw@x, round(raw@x))), "RNA_raw_counts integer-valued")
check(isTRUE(all.equal(cb@x, round(cb@x))), "RNA_cellbender_counts integer-valued")
d <- raw - cb
check(length(d@x) == 0 || min(d@x) >= 0, "RNA_cellbender_counts <= RNA_raw_counts element-wise")

# ─── Barcodes / features vs CellRanger ──────────────────────────────────────
cat("\n[cellranger]\n")
cr_barcodes <- readLines(gzfile(file.path(mtx_dir, "barcodes.tsv.gz")))
cr_barcodes <- cr_barcodes[nzchar(cr_barcodes)]
cr_feat_tbl <- read.delim(gzfile(file.path(mtx_dir, "features.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)
cr_features <- cr_feat_tbl[[2]]
cr_features_unique <- make.unique(cr_features)

check(identical(colnames(obj), cr_barcodes), sprintf("colnames match CellRanger barcodes (%d)", length(cr_barcodes)))
check(identical(rownames(obj), cr_features_unique), sprintf("rownames match CellRanger features (%d)", length(cr_features)))

cr_mtx   <- readMM(gzfile(file.path(mtx_dir, "matrix.mtx.gz")))
cr_total <- sum(cr_mtx)
raw_total <- sum(raw)
check(isTRUE(all.equal(cr_total, raw_total)), sprintf("total raw counts match CellRanger (%.0f vs %.0f)", raw_total, cr_total))

# ─── Souporcell ─────────────────────────────────────────────────────────────
cat("\n[souporcell]\n")
meta <- obj[[]]
for (col in c("souporcell_status", "souporcell_assignment")) {
    check(col %in% colnames(meta), sprintf("metadata column '%s' present", col))
}
status <- as.character(meta$souporcell_status)
assign <- as.character(meta$souporcell_assignment)
names(status) <- names(assign) <- colnames(obj)

if (souporcell_dir == "NONE") {
    check(all(is.na(status)) && all(is.na(assign)), "k=1: souporcell columns all NA")
} else {
    clusters <- file.path(souporcell_dir, "clusters.tsv")
    if (file.exists(clusters)) {
        sc_df <- read.delim(clusters, stringsAsFactors = FALSE, colClasses = "character")
        rownames(sc_df) <- sc_df$barcode
        shared <- intersect(colnames(obj), rownames(sc_df))
        check(length(shared) > 0, sprintf("souporcell barcodes overlap object (%d)", length(shared)))
        check(all(unique(na.omit(status)) %in% VALID_STATUS),
              sprintf("status values valid: %s", paste(sort(unique(na.omit(status))), collapse = ", ")))
        check(identical(unname(status[shared]), sc_df[shared, "status"]), "status matches clusters.tsv for shared barcodes")
        check(identical(unname(assign[shared]), sc_df[shared, "assignment"]), "assignment matches clusters.tsv for shared barcodes")
        absent <- setdiff(colnames(obj), rownames(sc_df))
        check(all(is.na(status[absent])), "barcodes absent from clusters.tsv are NA")
    } else {
        check(all(is.na(status)), "clusters.tsv missing: souporcell columns all NA")
    }
}

# ─── Sample-sheet metadata ──────────────────────────────────────────────────
cat("\n[metadata]\n")
sheet_meta <- fromJSON(metadata_json)
for (col in names(sheet_meta)) {
    present <- col %in% colnames(meta)
    check(present, sprintf("metadata column '%s' present", col))
    if (present) {
        vals <- unique(na.omit(as.character(meta[[col]])))
        expected <- sheet_meta[[col]]
        if (is.null(expected) || identical(expected, "")) {
            check(length(vals) == 0, sprintf("'%s' is NA (empty in sample sheet)", col))
        } else {
            check(length(vals) == 1 && vals[1] == as.character(expected), sprintf("'%s' == '%s' for all cells", col, expected))
        }
    }
}
check(as.character(meta$library_id[1]) == library_id, "library_id matches")

# ─── Summary ────────────────────────────────────────────────────────────────
status_counts <- table(ifelse(is.na(status), "NA", status))
assign_counts <- table(ifelse(is.na(assign), "NA", assign))
summary <- list(
    library_id                   = library_id,
    format                       = "seurat",
    n_cells                      = ncol(obj),
    n_features                   = nrow(obj),
    sum_raw_counts               = raw_total,
    sum_cellbender_counts        = sum(cb),
    nnz_raw                      = length(raw@x),
    nnz_cellbender               = length(cb@x),
    barcodes_md5                 = md5_of(colnames(obj)),
    features_md5                 = md5_of(cr_features),
    souporcell_status_counts     = as.list(setNames(as.integer(status_counts), names(status_counts))),
    souporcell_assignment_counts = as.list(setNames(as.integer(assign_counts), names(assign_counts))),
    sample_sheet                 = lapply(sheet_meta, function(v) if (is.null(v)) "" else as.character(v))
)
write_json(summary, paste0(library_id, "_seurat_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)

cat(sprintf("\n=== %d failure(s) ===\n", length(failures)))
for (f in failures) cat("  -", f, "\n")
quit(status = if (length(failures) > 0) 1 else 0)
