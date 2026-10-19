#!/usr/bin/env python3
"""
Validate an xRAFT-scRNA AnnData object against its inputs.

Checks (all must pass, exit 1 otherwise):
  structure  : X, layers RNA_raw_counts / RNA_cellbender_counts present, same shape
  identity   : X == RNA_raw_counts exactly
  values     : counts are finite, non-negative, integer-valued
  cellbender : 0 <= RNA_cellbender_counts <= RNA_raw_counts element-wise
  barcodes   : obs_names == CellRanger filtered barcodes.tsv.gz (same order)
  features   : var_names == CellRanger features (made unique, same order)
  counts     : sum(RNA_raw_counts) == sum(CellRanger matrix.mtx.gz)
  souporcell : status/assignment columns present; k>1 => assigned cells exist
               and status values are valid; k==1 => all NA
  metadata   : every sample-sheet column present in obs with exactly one value

Writes <library_id>_anndata_summary.json for cross-format comparison.
"""

import gzip
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io
import anndata as ad

VALID_STATUS = {"singlet", "doublet", "unassigned"}
failures = []


def check(cond, msg):
    tag = "PASS" if cond else "FAIL"
    print(f"  [{tag}] {msg}", flush=True)
    if not cond:
        failures.append(msg)


def md5_of(items):
    return hashlib.md5("\n".join(map(str, items)).encode()).hexdigest()


def make_unique(names):
    # scanpy-style: duplicate names get -1, -2, ...
    return list(ad.utils.make_index_unique(pd.Index(names)))


if len(sys.argv) != 6:
    sys.exit("Usage: validate_anndata.py <library_id> <h5ad> <cellranger_outs> <souporcell_dir|NONE> <metadata_json>")

library_id, h5ad_path, cellranger_outs, souporcell_dir, metadata_json = sys.argv[1:6]
mtx_dir = os.path.join(cellranger_outs, "filtered_feature_bc_matrix")

print(f"=== Validating AnnData: {library_id} ===")
print(f"h5ad            : {h5ad_path}")
print(f"CellRanger outs : {cellranger_outs}")
print(f"Souporcell dir  : {souporcell_dir}")

# ─── Load ──────────────────────────────────────────────────────────────────
adata = ad.read_h5ad(h5ad_path)
print(f"\nLoaded: {adata.n_obs} cells x {adata.n_vars} features")

# ─── Structure ─────────────────────────────────────────────────────────────
print("\n[structure]")
for layer in ("RNA_raw_counts", "RNA_cellbender_counts"):
    check(layer in adata.layers, f"layer '{layer}' present")
if failures:
    sys.exit("Missing layers; aborting further checks.")

X   = sp.csr_matrix(adata.X)
raw = sp.csr_matrix(adata.layers["RNA_raw_counts"])
cb  = sp.csr_matrix(adata.layers["RNA_cellbender_counts"])
check(X.shape == raw.shape == cb.shape == (adata.n_obs, adata.n_vars), "all layers share the obs x var shape")
check(adata.n_obs > 0 and adata.n_vars > 0, "non-empty object")

# ─── Identity / values ─────────────────────────────────────────────────────
print("\n[values]")
check((X != raw).nnz == 0, "X is identical to RNA_raw_counts")
check(np.isfinite(raw.data).all() and np.isfinite(cb.data).all(), "all counts finite")
check((raw.data >= 0).all() and (cb.data >= 0).all(), "all counts non-negative")
check(np.allclose(raw.data, np.round(raw.data)), "RNA_raw_counts integer-valued")
check(np.allclose(cb.data, np.round(cb.data)), "RNA_cellbender_counts integer-valued")
diff = raw - cb
check(diff.nnz == 0 or diff.data.min() >= 0, "RNA_cellbender_counts <= RNA_raw_counts element-wise")

# ─── Barcodes / features vs CellRanger ─────────────────────────────────────
print("\n[cellranger]")
with gzip.open(os.path.join(mtx_dir, "barcodes.tsv.gz"), "rt") as fh:
    cr_barcodes = [l.strip() for l in fh if l.strip()]
with gzip.open(os.path.join(mtx_dir, "features.tsv.gz"), "rt") as fh:
    cr_features = [l.rstrip("\n").split("\t")[1] for l in fh if l.strip()]
cr_features_unique = make_unique(cr_features)

check(list(adata.obs_names) == cr_barcodes, f"obs_names match CellRanger barcodes ({len(cr_barcodes)})")
check(list(adata.var_names) == cr_features_unique, f"var_names match CellRanger features ({len(cr_features)})")

cr_mtx = scipy.io.mmread(os.path.join(mtx_dir, "matrix.mtx.gz"))
cr_total = float(cr_mtx.sum())
raw_total = float(raw.sum())
check(np.isclose(cr_total, raw_total), f"total raw counts match CellRanger ({raw_total:.0f} vs {cr_total:.0f})")

# ─── Souporcell ────────────────────────────────────────────────────────────
print("\n[souporcell]")
for col in ("souporcell_status", "souporcell_assignment"):
    check(col in adata.obs.columns, f"obs column '{col}' present")

status = adata.obs["souporcell_status"].astype("string")
assign = adata.obs["souporcell_assignment"].astype("string")
status_counts = status.fillna("NA").value_counts().to_dict()

if souporcell_dir == "NONE":
    check(status.isna().all() and assign.isna().all(), "k=1: souporcell columns all NA")
else:
    clusters = os.path.join(souporcell_dir, "clusters.tsv")
    if os.path.isfile(clusters):
        sc_df = pd.read_csv(clusters, sep="\t", dtype=str).set_index("barcode")
        shared = adata.obs_names.intersection(sc_df.index)
        check(len(shared) > 0, f"souporcell barcodes overlap object ({len(shared)})")
        check(set(status.dropna().unique()) <= VALID_STATUS, f"status values valid: {sorted(status.dropna().unique())}")
        check((status.loc[shared].values == sc_df.loc[shared, "status"].values).all(), "status matches clusters.tsv for shared barcodes")
        check((assign.loc[shared].values == sc_df.loc[shared, "assignment"].values).all(), "assignment matches clusters.tsv for shared barcodes")
        check(status.loc[adata.obs_names.difference(sc_df.index)].isna().all(), "barcodes absent from clusters.tsv are NA")
    else:
        check(status.isna().all(), "clusters.tsv missing: souporcell columns all NA")

# ─── Sample-sheet metadata ─────────────────────────────────────────────────
print("\n[metadata]")
with open(metadata_json) as fh:
    sheet_meta = json.load(fh)
for col, value in sheet_meta.items():
    present = col in adata.obs.columns
    check(present, f"obs column '{col}' present")
    if present:
        vals = adata.obs[col].astype("string").dropna().unique()
        expected = None if value in (None, "") else str(value)
        if expected is None:
            check(len(vals) == 0, f"'{col}' is NA (empty in sample sheet)")
        else:
            check(len(vals) == 1 and vals[0] == expected, f"'{col}' == '{expected}' for all cells")
check(str(adata.obs["library_id"].iloc[0]) == library_id, "library_id matches")

# ─── Summary ───────────────────────────────────────────────────────────────
summary = {
    "library_id": library_id,
    "format": "anndata",
    "n_cells": int(adata.n_obs),
    "n_features": int(adata.n_vars),
    "sum_raw_counts": raw_total,
    "sum_cellbender_counts": float(cb.sum()),
    "nnz_raw": int(raw.nnz),
    "nnz_cellbender": int(cb.nnz),
    "barcodes_md5": md5_of(adata.obs_names),
    "features_md5": md5_of(cr_features),
    "souporcell_status_counts": {k: int(v) for k, v in status_counts.items()},
    "souporcell_assignment_counts": {str(k): int(v) for k, v in assign.fillna("NA").value_counts().to_dict().items()},
    "sample_sheet": {k: ("" if v is None else str(v)) for k, v in sheet_meta.items()},
}
with open(f"{library_id}_anndata_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2, sort_keys=True)

print(f"\n=== {len(failures)} failure(s) ===")
for f in failures:
    print(f"  - {f}")
sys.exit(1 if failures else 0)
