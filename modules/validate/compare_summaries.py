#!/usr/bin/env python3
"""
Compare the Seurat and AnnData validation summaries for one library.
Standard library only, so it runs on the bare cluster python3.

Exit 1 if any field disagrees.
"""

import json
import math
import sys

if len(sys.argv) != 4:
    sys.exit("Usage: compare_summaries.py <library_id> <seurat_summary.json> <anndata_summary.json>")

library_id, seurat_json, anndata_json = sys.argv[1:4]

with open(seurat_json) as fh:
    s = json.load(fh)
with open(anndata_json) as fh:
    a = json.load(fh)

failures = []


def check(cond, msg):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        failures.append(msg)


def same_num(x, y):
    return math.isclose(float(x), float(y), rel_tol=1e-9, abs_tol=0.5)


print(f"=== Comparing Seurat vs AnnData: {library_id} ===")
check(s["library_id"] == a["library_id"] == library_id, "library_id")
for key in ("n_cells", "n_features", "nnz_raw", "nnz_cellbender"):
    check(s[key] == a[key], f"{key}: {s[key]} vs {a[key]}")
for key in ("sum_raw_counts", "sum_cellbender_counts"):
    check(same_num(s[key], a[key]), f"{key}: {s[key]} vs {a[key]}")
for key in ("barcodes_md5", "features_md5"):
    check(s[key] == a[key], f"{key} identical")
for key in ("souporcell_status_counts", "souporcell_assignment_counts", "sample_sheet"):
    check(s[key] == a[key], f"{key}: {s[key]} vs {a[key]}")

print(f"\n=== {len(failures)} disagreement(s) ===")
for f in failures:
    print(f"  - {f}")
sys.exit(1 if failures else 0)
