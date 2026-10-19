#!/bin/bash
# ─── xRAFT-scRNA unit tests ───────────────────────────────────────────────────
# Runs the object-creation scripts and their validators on synthetic data,
# for k>1 (Souporcell present) and k=1 (no Souporcell), then checks that the
# Seurat and AnnData objects agree.
#
# Usage:
#   tests/run_tests.sh                 # use the pipeline containers (apptainer)
#   tests/run_tests.sh --local         # use python3 / Rscript on PATH
#   tests/run_tests.sh --python-only   # skip the R side (no Seurat available)
#
# Containers are pulled into ./containers if missing (same SIF names the
# pipeline uses, so a later pipeline run reuses them).

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${HERE}/.." && pwd)"
OUT="${HERE}/output"
MODE="container"
RUN_R=true

for arg in "$@"; do
    case "$arg" in
        --local)       MODE="local" ;;
        --python-only) RUN_R=false ;;
        *) echo "Unknown option: $arg"; exit 2 ;;
    esac
done

# ── Interpreters ────────────────────────────────────────────────────────────
# Images come from nextflow.config so tests run in exactly the pinned containers.
cfg_image() { grep -E "^\s*$1\s*=" "${ROOT}/nextflow.config" | sed -E 's/.*"([^"]+)".*/\1/'; }
sif_name()  { local img="$2"; local tag="${img##*:}"; echo "${ROOT}/containers/$1_${tag//[^A-Za-z0-9._-]/_}.sif"; }

if [ "$MODE" = "container" ]; then
    SCANPY_IMG=$(cfg_image scanpy_container); SCANPY_SIF=$(sif_name scanpy "$SCANPY_IMG")
    SEURAT_IMG=$(cfg_image seurat_container); SEURAT_SIF=$(sif_name seurat "$SEURAT_IMG")
    mkdir -p "${ROOT}/containers"
    export APPTAINER_CACHEDIR="${OUT}/.apptainer_cache" APPTAINER_TMPDIR="${OUT}/.apptainer_tmp"
    mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
    if [ ! -f "$SCANPY_SIF" ]; then
        apptainer pull "$SCANPY_SIF" "$SCANPY_IMG"
        bash "${ROOT}/modules/containers/record_digest.sh" "$SCANPY_SIF" "$SCANPY_IMG" "$APPTAINER_CACHEDIR"
    fi
    PY="apptainer exec --bind ${ROOT} ${SCANPY_SIF} python3"
    if $RUN_R; then
        if [ ! -f "$SEURAT_SIF" ]; then
            apptainer pull "$SEURAT_SIF" "$SEURAT_IMG"
            bash "${ROOT}/modules/containers/record_digest.sh" "$SEURAT_SIF" "$SEURAT_IMG" "$APPTAINER_CACHEDIR"
        fi
        R="apptainer exec --bind ${ROOT} ${SEURAT_SIF} Rscript"
    fi
else
    PY="python3"
    R="Rscript"
fi

rm -rf "$OUT"
mkdir -p "$OUT"
FAILED=0

run() {  # run <label> <command...>
    local label="$1"; shift
    echo ""
    echo "────────────────────────────────────────────────────────────"
    echo "  $label"
    echo "────────────────────────────────────────────────────────────"
    if "$@"; then
        echo "  ✔ $label"
    else
        echo "  ✘ $label"
        FAILED=$((FAILED + 1))
    fi
}

for K in 2 1; do
    ID="TEST_K${K}"
    DATA="${OUT}/data_k${K}"
    WORK="${OUT}/work_k${K}"
    mkdir -p "$WORK"

    echo ""
    echo "============================================================"
    echo "  Synthetic dataset: ${ID} (souporcell_k = ${K})"
    echo "============================================================"
    SPEC=$($PY "${HERE}/make_synthetic_data.py" --out "$DATA" --id "$ID" --k "$K")
    echo "$SPEC"
    CR_OUTS=$(echo "$SPEC" | $PY -c 'import json,sys; print(json.load(sys.stdin)["cellranger_outs"])')
    CB_DIR=$(echo "$SPEC"  | $PY -c 'import json,sys; print(json.load(sys.stdin)["cellbender_dir"])')
    SC_DIR=$(echo "$SPEC"  | $PY -c 'import json,sys; print(json.load(sys.stdin)["souporcell_dir"])')
    META=$(echo "$SPEC"    | $PY -c 'import json,sys; print(json.load(sys.stdin)["metadata_json"])')

    # ── AnnData ─────────────────────────────────────────────────────────
    run "create_anndata.py (${ID})" bash -c "cd '$WORK' && $PY '${ROOT}/modules/create_anndata_object/create_anndata.py' '$ID' '$CR_OUTS' '$CB_DIR' '$SC_DIR' '$META' > create_anndata.log"
    run "validate_anndata.py (${ID})" bash -c "cd '$WORK' && $PY '${ROOT}/modules/validate/validate_anndata.py' '$ID' '${ID}.h5ad' '$CR_OUTS' '$SC_DIR' '$META' | tee validate_anndata.log"

    # ── Seurat ──────────────────────────────────────────────────────────
    if $RUN_R; then
        run "create_seurat.R (${ID})" bash -c "cd '$WORK' && $R '${ROOT}/modules/create_seurat_object/create_seurat.R' '$ID' '$CR_OUTS' '$CB_DIR' '$SC_DIR' '$META' > create_seurat.log"
        run "validate_seurat.R (${ID})" bash -c "cd '$WORK' && $R '${ROOT}/modules/validate/validate_seurat.R' '$ID' '${ID}_seurat.rds' '$CR_OUTS' '$SC_DIR' '$META' | tee validate_seurat.log"
        run "compare_summaries.py (${ID})" bash -c "cd '$WORK' && python3 '${ROOT}/modules/validate/compare_summaries.py' '$ID' '${ID}_seurat_summary.json' '${ID}_anndata_summary.json' | tee compare.log"
    fi
done

echo ""
echo "============================================================"
if [ "$FAILED" -eq 0 ]; then
    echo "  ALL TESTS PASSED"
else
    echo "  ${FAILED} TEST(S) FAILED — see ${OUT}/work_k*/ logs"
fi
echo "============================================================"
exit $FAILED
