#!/usr/bin/env Rscript

# ─── Create Seurat Object ───────────────────────────────────────────────────
# Builds a Seurat object with:
#   - RNA assay (Assay5) with layers:
#       counts / RNA_raw_counts : CellRanger filtered counts
#       RNA_cellbender_counts   : CellBender corrected counts (aligned to CellRanger grid)
#   - Souporcell metadata  : souporcell_status and souporcell_assignment from clusters.tsv
#   - Sample-sheet metadata: every non-path column of the sample sheet, one value per cell

# Runs inside the satijalab/seurat container (R 4.5, Seurat 5, hdf5r, jsonlite).
cat("R", R.version.string, "| Seurat", as.character(packageVersion("Seurat")), "\n\n")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
    stop("Usage: create_seurat.R <library_id> <cellranger_outs> <cellbender_dir> <souporcell_dir> <metadata_json>")
}

library_id      <- args[1]
cellranger_outs <- args[2]
cellbender_dir  <- args[3]
souporcell_dir  <- args[4]
metadata_json   <- args[5]

# Auto-discover CellBender filtered h5 file
cb_h5_candidates <- list.files(cellbender_dir, pattern = "filtered.*\\.h5$", full.names = TRUE)
if (length(cb_h5_candidates) == 0) {
    cb_h5_candidates <- list.files(cellbender_dir, pattern = "\\.h5$", full.names = TRUE)
}
if (length(cb_h5_candidates) == 0) {
    stop("ERROR: No .h5 files found in CellBender directory: ", cellbender_dir,
         "\n  Contents: ", paste(list.files(cellbender_dir), collapse = ", "))
}
cellbender_h5 <- cb_h5_candidates[1]

cat("=== Creating Seurat object for:", library_id, "===\n")
cat("CellRanger outs :", cellranger_outs, "\n")
cat("CellBender dir  :", cellbender_dir, "\n")
cat("CellBender h5   :", cellbender_h5, "\n")
cat("Souporcell dir  :", souporcell_dir, "\n")
cat("Metadata JSON   :", metadata_json, "\n")

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(hdf5r)
    library(jsonlite)
})

# ─── CellBender h5 reader (10x v3 layout; replaces scCustomize::Read_CellBender_h5_Mat)
read_cellbender_h5 <- function(file_name) {
    h5 <- H5File$new(file_name, mode = "r")
    on.exit(h5$close_all())
    g <- h5[["matrix"]]
    feature_names <- g[["features/name"]][]
    feature_ids   <- g[["features/id"]][]
    barcodes      <- g[["barcodes"]][]
    dims          <- g[["shape"]][]
    mat <- sparseMatrix(
        i        = g[["indices"]][] + 1,
        p        = g[["indptr"]][],
        x        = as.numeric(g[["data"]][]),
        dims     = dims,
        dimnames = list(make.unique(feature_names), barcodes),
        repr     = "C"
    )
    attr(mat, "feature_ids") <- feature_ids
    mat
}

# ─── 1. RNA assay (Assay5) from CellRanger filtered matrix ─────────────────
# Layer layout mirrors the AnnData object:
#   counts               : CellRanger raw counts (Seurat's required default layer, = AnnData X)
#   RNA_raw_counts       : CellRanger raw counts
#   RNA_cellbender_counts: CellBender corrected counts, aligned to the full feature x barcode grid
cat("\n[1/5] Loading CellRanger filtered matrix...\n")
cr_matrix <- Read10X(data.dir = file.path(cellranger_outs, "filtered_feature_bc_matrix"))
seurat_obj <- CreateSeuratObject(
    counts   = cr_matrix,
    project  = library_id,
    assay    = "RNA"
)
cat("  Barcodes:", ncol(seurat_obj), "| Features:", nrow(seurat_obj), "\n")

cr_barcodes <- colnames(seurat_obj)
cr_features <- rownames(seurat_obj)

LayerData(seurat_obj, assay = "RNA", layer = "RNA_raw_counts") <- LayerData(seurat_obj, assay = "RNA", layer = "counts")
cat("  RNA_raw_counts layer added.\n")

# ─── 2. RNA_cellbender_counts layer from CellBender output ──────────────────
cat("\n[2/5] Loading CellBender output...\n")
cb_matrix <- read_cellbender_h5(cellbender_h5)

shared_barcodes <- intersect(cr_barcodes, colnames(cb_matrix))
shared_features <- intersect(cr_features, rownames(cb_matrix))
cat("  CellBender total barcodes:", ncol(cb_matrix), "\n")
cat("  Shared with CellRanger   :", length(shared_barcodes), "\n")
cat("  Shared features          :", length(shared_features), "\n")

if (length(shared_barcodes) == 0) {
    stop("ERROR: No overlapping barcodes between CellRanger and CellBender outputs.")
}

missing_barcodes <- setdiff(cr_barcodes, shared_barcodes)
missing_features <- setdiff(cr_features, shared_features)
if (length(missing_barcodes) > 0)
    cat("  Barcodes in CellRanger but not CellBender:", length(missing_barcodes), "(filled with zeros)\n")
if (length(missing_features) > 0)
    cat("  Features in CellRanger but not CellBender:", length(missing_features), "(filled with zeros)\n")

# Build a full feature x barcode sparse matrix in CellRanger order, zero-filled
cb_sub <- as(cb_matrix[shared_features, shared_barcodes, drop = FALSE], "TsparseMatrix")
cb_aligned <- sparseMatrix(
    i        = match(rownames(cb_sub), cr_features)[cb_sub@i + 1],
    j        = match(colnames(cb_sub), cr_barcodes)[cb_sub@j + 1],
    x        = cb_sub@x,
    dims     = c(length(cr_features), length(cr_barcodes)),
    dimnames = list(cr_features, cr_barcodes),
    repr     = "C"
)

LayerData(seurat_obj, assay = "RNA", layer = "RNA_cellbender_counts") <- cb_aligned
cat("  RNA_cellbender_counts layer added.\n")

# ─── 3. Souporcell metadata ─────────────────────────────────────────────────
cat("\n[3/5] Loading Souporcell clusters.tsv...\n")

souporcell_status     <- rep(NA_character_, length(cr_barcodes))
souporcell_assignment <- rep(NA_character_, length(cr_barcodes))
names(souporcell_status)     <- cr_barcodes
names(souporcell_assignment) <- cr_barcodes

if (souporcell_dir == "NONE") {
    cat("  Souporcell not run for this sample (k=1), setting metadata to NA.\n")
} else {
    clusters_file <- file.path(souporcell_dir, "clusters.tsv")

    if (file.exists(clusters_file)) {
        sc_clusters <- read.delim(clusters_file, stringsAsFactors = FALSE)
        cat("  Souporcell barcodes:", nrow(sc_clusters), "\n")

        rownames(sc_clusters) <- sc_clusters$barcode
        shared_sc <- intersect(cr_barcodes, rownames(sc_clusters))
        cat("  Shared with CellRanger:", length(shared_sc), "\n")

        souporcell_status[shared_sc]     <- as.character(sc_clusters[shared_sc, "status"])
        souporcell_assignment[shared_sc] <- as.character(sc_clusters[shared_sc, "assignment"])
        cat("  Souporcell metadata added (souporcell_status, souporcell_assignment).\n")
    } else {
        warning("Souporcell clusters.tsv not found at: ", clusters_file)
    }
}

seurat_obj$souporcell_status     <- souporcell_status
seurat_obj$souporcell_assignment <- souporcell_assignment

# ─── 4. Sample-sheet metadata ───────────────────────────────────────────────
cat("\n[4/5] Adding sample-sheet metadata...\n")
sheet_meta <- fromJSON(metadata_json)
for (col in names(sheet_meta)) {
    value <- sheet_meta[[col]]
    if (is.null(value) || identical(value, "")) value <- NA_character_
    seurat_obj[[col]] <- as.character(value)
    cat("  ", col, ":", as.character(value), "\n")
}

# ─── 5. Save Seurat object ──────────────────────────────────────────────────
cat("\n[5/5] Saving Seurat object...\n")
output_file <- paste0(library_id, "_seurat.rds")
saveRDS(seurat_obj, file = output_file)
cat("  Saved to:", output_file, "\n")

cat("\n=== Summary ===\n")
cat("Assay          :", DefaultAssay(seurat_obj), "(", class(seurat_obj[["RNA"]])[1], ")\n")
cat("Layers         :", paste(Layers(seurat_obj[["RNA"]]), collapse = ", "), "\n")
cat("Cells          :", ncol(seurat_obj), "\n")
cat("Features       :", nrow(seurat_obj[["RNA"]]), "\n")
cat("Metadata cols  :", paste(colnames(seurat_obj[[]]), collapse = ", "), "\n")
cat("Done.\n")
