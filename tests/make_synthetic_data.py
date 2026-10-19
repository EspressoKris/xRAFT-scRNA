#!/usr/bin/env python3
"""
Build a small synthetic input set that mimics what the pipeline hands to the
object-creation scripts, so create_seurat.R / create_anndata.py and their
validators can be exercised in seconds without Cell Ranger, CellBender or
Souporcell.

Layout written under <out_dir>:
  cellranger/<id>/outs/filtered_feature_bc_matrix/{matrix.mtx.gz,barcodes.tsv.gz,features.tsv.gz}
  cellbender/<id>/<id>_cellbender_output_filtered.h5      (10x v3-layout h5)
  souporcell/<id>/clusters.tsv                           (only when k > 1)
  <id>_metadata.json

Deliberate edge cases:
  - a slice of CellRanger barcodes is absent from CellBender   (zero-fill path)
  - CellBender carries extra barcodes not in CellRanger         (must be dropped)
  - one duplicated gene symbol                                  (make-unique path)
  - CellBender counts are a random subset of raw counts         (<= raw holds)
  - some CellRanger barcodes are missing from clusters.tsv      (NA path)
  - one metadata column is empty                                (NA metadata)
"""

import argparse
import gzip
import json
import os

import h5py
import numpy as np
import scipy.io
import scipy.sparse as sp


def write_10x_mtx(out_dir, mat, barcodes, features):
    os.makedirs(out_dir, exist_ok=True)
    with gzip.open(os.path.join(out_dir, "matrix.mtx.gz"), "wb") as fh:
        scipy.io.mmwrite(fh, sp.coo_matrix(mat), field="integer")
    with gzip.open(os.path.join(out_dir, "barcodes.tsv.gz"), "wt") as fh:
        fh.write("\n".join(barcodes) + "\n")
    with gzip.open(os.path.join(out_dir, "features.tsv.gz"), "wt") as fh:
        for fid, name in features:
            fh.write(f"{fid}\t{name}\tGene Expression\n")


def write_cellbender_h5(path, mat, barcodes, features):
    """CellBender 'filtered' output: 10x v3 layout under /matrix plus extras."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    csc = sp.csc_matrix(mat)
    with h5py.File(path, "w") as h5:
        g = h5.create_group("matrix")
        g.create_dataset("data", data=csc.data.astype(np.float32))
        g.create_dataset("indices", data=csc.indices.astype(np.int64))
        g.create_dataset("indptr", data=csc.indptr.astype(np.int64))
        g.create_dataset("shape", data=np.array(csc.shape, dtype=np.int32))
        g.create_dataset("barcodes", data=np.array(barcodes, dtype="S"))
        f = g.create_group("features")
        f.create_dataset("id", data=np.array([fid for fid, _ in features], dtype="S"))
        f.create_dataset("name", data=np.array([name for _, name in features], dtype="S"))
        f.create_dataset("feature_type", data=np.array(["Gene Expression"] * len(features), dtype="S"))
        f.create_dataset("genome", data=np.array(["GRCh38"] * len(features), dtype="S"))
        # CellBender-specific extras (ignored by readers, present in real files)
        g.create_dataset("latent_cell_probability", data=np.random.rand(len(barcodes)).astype(np.float32))
        h5.attrs["cellbender_version"] = "synthetic"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--id", default="TEST_1")
    ap.add_argument("--k", type=int, default=2, help="souporcell_k (1 = no souporcell)")
    ap.add_argument("--cells", type=int, default=300)
    ap.add_argument("--genes", type=int, default=120)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    out = args.out
    lib = args.id

    # ── features (one duplicated symbol on purpose) ─────────────────────
    features = [(f"ENSG{i:08d}", f"GENE{i}") for i in range(args.genes)]
    features[5] = ("ENSG_DUP", "GENE3")  # duplicate of GENE3

    # ── CellRanger filtered matrix ──────────────────────────────────────
    cr_barcodes = [f"{''.join(rng.choice(list('ACGT'), 16))}-1" for _ in range(args.cells)]
    cr_barcodes = list(dict.fromkeys(cr_barcodes))  # ensure unique
    raw = rng.poisson(0.8, size=(args.genes, len(cr_barcodes)))
    raw[:, : len(cr_barcodes) // 3] += rng.poisson(2, size=(args.genes, len(cr_barcodes) // 3))
    write_10x_mtx(os.path.join(out, "cellranger", lib, "outs", "filtered_feature_bc_matrix"), raw, cr_barcodes, features)

    # ── CellBender: subset of raw, drop some CR barcodes, add extras ────
    keep = np.ones(len(cr_barcodes), dtype=bool)
    keep[-10:] = False                                   # 10 CR barcodes missing from CB
    cb_barcodes = [b for b, k in zip(cr_barcodes, keep) if k]
    cb = raw[:, keep].copy()
    removal = rng.random(cb.shape) < 0.3                  # ambient removal
    cb = np.where(removal, np.maximum(cb - 1, 0), cb)
    extra_bcs = [f"{''.join(rng.choice(list('ACGT'), 16))}-1" for _ in range(5)]
    cb = np.hstack([cb, rng.poisson(0.3, size=(args.genes, 5))])
    cb_barcodes = cb_barcodes + extra_bcs
    perm = rng.permutation(len(cb_barcodes))              # CellBender order differs from CR
    write_cellbender_h5(
        os.path.join(out, "cellbender", lib, f"{lib}_cellbender_output_filtered.h5"),
        cb[:, perm], [cb_barcodes[i] for i in perm], features,
    )

    # ── Souporcell ──────────────────────────────────────────────────────
    sc_dir = "NONE"
    if args.k > 1:
        sc_dir = os.path.join(out, "souporcell", lib)
        os.makedirs(sc_dir, exist_ok=True)
        n = len(cr_barcodes) - 7                          # 7 CR barcodes absent from clusters.tsv
        status = rng.choice(["singlet", "doublet", "unassigned"], size=n, p=[0.8, 0.15, 0.05])
        assignment = [
            (str(rng.integers(args.k)) if s == "singlet" else f"{rng.integers(args.k)}/{rng.integers(args.k)}")
            for s in status
        ]
        with open(os.path.join(sc_dir, "clusters.tsv"), "w") as fh:
            fh.write("barcode\tstatus\tassignment\tlog_prob_singleton\tlog_prob_doublet\n")
            for b, s, a in zip(cr_barcodes[:n], status, assignment):
                fh.write(f"{b}\t{s}\t{a}\t{-rng.random():.3f}\t{-rng.random():.3f}\n")

    # ── Sample-sheet metadata ───────────────────────────────────────────
    meta = {"library_id": lib, "experiment": "MF", "batch_10x": "202602", "souporcell_k": str(args.k), "notes": ""}
    with open(os.path.join(out, f"{lib}_metadata.json"), "w") as fh:
        json.dump(meta, fh)

    print(json.dumps({
        "library_id": lib,
        "cellranger_outs": os.path.join(out, "cellranger", lib, "outs"),
        "cellbender_dir": os.path.join(out, "cellbender", lib),
        "souporcell_dir": sc_dir,
        "metadata_json": os.path.join(out, f"{lib}_metadata.json"),
        "n_cells": len(cr_barcodes),
        "n_features": args.genes,
        "sum_raw_counts": int(raw.sum()),
    }))


if __name__ == "__main__":
    main()
