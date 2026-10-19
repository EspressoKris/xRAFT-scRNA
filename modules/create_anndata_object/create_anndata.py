#!/usr/bin/env python3
"""
Create AnnData object — Python counterpart of create_seurat.R.

Builds an AnnData object with:
  - X / layers["RNA_raw_counts"] : CellRanger filtered counts (raw UMIs)
  - layers["RNA_cellbender_counts"]     : CellBender corrected counts (subset to CellRanger barcodes)
  - obs souporcell_status /
    souporcell_assignment     : from Souporcell clusters.tsv (NA when not run)
  - obs <sample-sheet column> : every non-path column of the sample sheet
  - uns["sample_sheet"]       : the same metadata as a dict
"""

import json
import os
import sys
import glob

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc


def log(msg=""):
    print(msg, flush=True)


# ─── Arguments ──────────────────────────────────────────────────────────────
if len(sys.argv) != 6:
    sys.exit("Usage: create_anndata.py <library_id> <cellranger_outs> <cellbender_dir> <souporcell_dir> <metadata_json>")

library_id, cellranger_outs, cellbender_dir, souporcell_dir, metadata_json = sys.argv[1:6]

# Auto-discover CellBender filtered h5 file
cb_h5_candidates = sorted(glob.glob(os.path.join(cellbender_dir, "*filtered*.h5")))
if not cb_h5_candidates:
    cb_h5_candidates = sorted(glob.glob(os.path.join(cellbender_dir, "*.h5")))
if not cb_h5_candidates:
    sys.exit(f"ERROR: No .h5 files found in CellBender directory: {cellbender_dir}\n"
             f"  Contents: {', '.join(os.listdir(cellbender_dir))}")
cellbender_h5 = cb_h5_candidates[0]

log(f"=== Creating AnnData object for: {library_id} ===")
log(f"CellRanger outs : {cellranger_outs}")
log(f"CellBender dir  : {cellbender_dir}")
log(f"CellBender h5   : {cellbender_h5}")
log(f"Souporcell dir  : {souporcell_dir}")
log(f"Metadata JSON   : {metadata_json}")
log(f"scanpy {sc.__version__} | anndata {ad.__version__}")

# ─── 1. Counts from CellRanger filtered matrix ──────────────────────────────
log("\n[1/5] Loading CellRanger filtered matrix...")
adata = sc.read_10x_mtx(
    os.path.join(cellranger_outs, "filtered_feature_bc_matrix"),
    var_names="gene_symbols",
    cache=False,
)
adata.var_names_make_unique()
adata.X = sp.csr_matrix(adata.X)
adata.layers["RNA_raw_counts"] = adata.X.copy()
adata.obs["library_id"] = library_id
log(f"  Barcodes: {adata.n_obs} | Features: {adata.n_vars}")

cr_barcodes = adata.obs_names

# ─── 2. CellBender layer ────────────────────────────────────────────────────
log("\n[2/5] Loading CellBender output...")
# CellBender writes a 10x-style h5 (plus extra fields); scanpy's reader handles it
cb = sc.read_10x_h5(cellbender_h5, gex_only=False)
cb.var_names_make_unique()
cb.X = sp.csr_matrix(cb.X)

shared_barcodes = cr_barcodes.intersection(cb.obs_names)
log(f"  CellBender total barcodes: {cb.n_obs}")
log(f"  Shared with CellRanger   : {len(shared_barcodes)}")

if len(shared_barcodes) == 0:
    sys.exit("ERROR: No overlapping barcodes between CellRanger and CellBender outputs.")

shared_features = adata.var_names.intersection(cb.var_names)
log(f"  Shared features          : {len(shared_features)}")

missing_barcodes = cr_barcodes.difference(cb.obs_names)
if len(missing_barcodes) > 0:
    log(f"  Barcodes in CellRanger but not CellBender: {len(missing_barcodes)} (filled with zeros)")

# Build a layer aligned to adata (obs x var); zero-fill anything not in CellBender
cb_sub = sp.coo_matrix(cb[shared_barcodes, shared_features].X)
obs_idx = adata.obs_names.get_indexer(shared_barcodes)
var_idx = adata.var_names.get_indexer(shared_features)

layer = sp.csr_matrix(
    (cb_sub.data.astype(np.float32), (obs_idx[cb_sub.row], var_idx[cb_sub.col])),
    shape=(adata.n_obs, adata.n_vars),
)
adata.layers["RNA_cellbender_counts"] = layer
log("  RNA_cellbender_counts layer added.")

# ─── 3. Souporcell metadata ─────────────────────────────────────────────────
log("\n[3/5] Loading Souporcell clusters.tsv...")
souporcell_status = pd.Series(pd.NA, index=cr_barcodes, dtype="string")
souporcell_assignment = pd.Series(pd.NA, index=cr_barcodes, dtype="string")

if souporcell_dir == "NONE":
    log("  Souporcell not run for this sample (k=1), setting metadata to NA.")
else:
    clusters_file = os.path.join(souporcell_dir, "clusters.tsv")
    if os.path.isfile(clusters_file):
        sc_clusters = pd.read_csv(clusters_file, sep="\t", dtype=str).set_index("barcode")
        log(f"  Souporcell barcodes: {len(sc_clusters)}")
        shared_sc = cr_barcodes.intersection(sc_clusters.index)
        log(f"  Shared with CellRanger: {len(shared_sc)}")
        souporcell_status.loc[shared_sc] = sc_clusters.loc[shared_sc, "status"].values
        souporcell_assignment.loc[shared_sc] = sc_clusters.loc[shared_sc, "assignment"].values
        log("  Souporcell metadata added (souporcell_status, souporcell_assignment).")
    else:
        log(f"  WARNING: Souporcell clusters.tsv not found at: {clusters_file}")

adata.obs["souporcell_status"] = souporcell_status.astype("category")
adata.obs["souporcell_assignment"] = souporcell_assignment.astype("category")

# ─── 4. Sample-sheet metadata ───────────────────────────────────────────────
log("\n[4/5] Adding sample-sheet metadata...")
with open(metadata_json) as fh:
    sheet_meta = json.load(fh)

for col, value in sheet_meta.items():
    value = str(value) if value not in (None, "") else pd.NA
    adata.obs[col] = pd.Series(value, index=cr_barcodes, dtype="string").astype("category")
    log(f"   {col}: {value}")

adata.uns["sample_sheet"] = {k: ("" if v is None else str(v)) for k, v in sheet_meta.items()}
adata.uns["pipeline"] = {
    "name": "xRAFT-scRNA",
    "cellranger_outs": cellranger_outs,
    "cellbender_h5": cellbender_h5,
    "souporcell_dir": souporcell_dir,
}

# ─── 5. Save ────────────────────────────────────────────────────────────────
log("\n[5/5] Saving AnnData object...")
output_file = f"{library_id}.h5ad"
adata.write_h5ad(output_file, compression="gzip")
log(f"  Saved to: {output_file}")

log("\n=== Summary ===")
log(f"Layers         : {', '.join(k for k in adata.layers.keys() if k is not None)}")
log(f"Cells          : {adata.n_obs}")
log(f"Features       : {adata.n_vars}")
log(f"obs columns    : {', '.join(adata.obs.columns)}")
log("Done.")
