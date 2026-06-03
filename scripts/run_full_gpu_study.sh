#!/usr/bin/env bash
# run_full_gpu_study.sh -- full GPU production pipeline for the 2D ABF-FR study.
#
# Safe by default: it PREPARES everything (reference, validation, benchmark,
# shards) and prints the plan, but only LAUNCHES the production shards when you
# pass --yes-production.  Validation failure aborts the launch (to protect GPU
# hours) unless you pass --skip-validation.
#
# Usage:
#   # Single A100
#   bash scripts/run_full_gpu_study.sh \
#       --config configs/two_dim_xi_x_production_gpu.yaml \
#       --stage production_gpu --num-gpus 1 --yes-production
#
#   # 4 or 8 GPUs
#   bash scripts/run_full_gpu_study.sh --config ... --stage production_gpu \
#       --num-gpus 4 --yes-production
#
#   # Dry plan only (no launch): omit --yes-production
#   bash scripts/run_full_gpu_study.sh --config ... --stage production_gpu --num-gpus 4
#
# Respects CUDA_VISIBLE_DEVICES (its GPU list is used; --num-gpus is capped to it).
# Per-shard logs: results/two_dim_xi_x/<stage>/logs/shard_NNN.log

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

CONFIG="configs/two_dim_xi_x_production_gpu.yaml"
STAGE="production_gpu"
NUM_GPUS=1
YES_PRODUCTION=""
SKIP_VALIDATION=""
SKIP_BENCHMARK=""
BENCH_STEPS=3000
PY="${PYTHON:-python3}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)          CONFIG="$2"; shift 2 ;;
        --stage)           STAGE="$2"; shift 2 ;;
        --num-gpus)        NUM_GPUS="$2"; shift 2 ;;
        --yes-production)  YES_PRODUCTION=1; shift ;;
        --skip-validation) SKIP_VALIDATION=1; shift ;;
        --skip-benchmark)  SKIP_BENCHMARK=1; shift ;;
        --bench-steps)     BENCH_STEPS="$2"; shift 2 ;;
        -h|--help)         sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

OUTPUT_ROOT="$($PY -c "import sys; sys.path.insert(0,'src'); from abffr import io_utils; print(io_utils.output_root(io_utils.load_config('$CONFIG')))")"
STAGE_DIR="$OUTPUT_ROOT/$STAGE"
LOG_DIR="$STAGE_DIR/logs"
SHARD_DIR="$STAGE_DIR/shards"
mkdir -p "$LOG_DIR" "$SHARD_DIR"

echo "========================================================"
echo "  ABF-FR 2D GPU study -- $STAGE"
echo "  config : $CONFIG"
echo "  gpus   : $NUM_GPUS (CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset})"
echo "  started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "========================================================"

# --------------------------------------------------------------------------- #
# Step 1 -- Reference (run if missing)
# --------------------------------------------------------------------------- #
REF_CSV="$OUTPUT_ROOT/reference/reference_profile.csv"
if [[ -f "$REF_CSV" ]]; then
    echo "-- Step 1: reference exists ($REF_CSV)"
else
    echo "-- Step 1: reference missing; computing..."
    $PY scripts/run_reference_2d.py --config "$CONFIG" || { echo "Reference FAILED"; exit 1; }
fi

# --------------------------------------------------------------------------- #
# Step 2 -- Validation (abort on failure unless skipped)
# --------------------------------------------------------------------------- #
if [[ -z "$SKIP_VALIDATION" ]]; then
    echo "-- Step 2: validation..."
    if ! $PY scripts/validate_torch_backend.py --config "$CONFIG"; then
        echo "!! VALIDATION FAILED -- refusing to launch production."
        echo "   Inspect $OUTPUT_ROOT/validation/validation_report.json,"
        echo "   fix the issue, or re-run with --skip-validation to override."
        exit 2
    fi
else
    echo "-- Step 2: validation SKIPPED (--skip-validation)"
fi

# --------------------------------------------------------------------------- #
# Step 3 -- Benchmark (informational)
# --------------------------------------------------------------------------- #
if [[ -z "$SKIP_BENCHMARK" ]]; then
    echo "-- Step 3: benchmark (n_steps=$BENCH_STEPS)..."
    $PY scripts/benchmark_backends.py --config "$CONFIG" --n-steps "$BENCH_STEPS" \
        --skip-numpy || echo "   (benchmark failed; continuing)"
else
    echo "-- Step 3: benchmark SKIPPED (--skip-benchmark)"
fi

# --------------------------------------------------------------------------- #
# Step 4 -- Resolve GPU list and create shards
# --------------------------------------------------------------------------- #
if [[ -n "${CUDA_VISIBLE_DEVICES:-}" ]]; then
    IFS=',' read -r -a GPU_LIST <<< "$CUDA_VISIBLE_DEVICES"
else
    GPU_LIST=(); for ((i=0;i<NUM_GPUS;i++)); do GPU_LIST+=("$i"); done
fi
N_GPU=${#GPU_LIST[@]}
(( NUM_GPUS < N_GPU )) && N_GPU=$NUM_GPUS
echo "-- Step 4: creating $N_GPU shard(s)..."
$PY scripts/make_gpu_shards.py --config "$CONFIG" --stage "$STAGE" --num-shards "$N_GPU"

# --------------------------------------------------------------------------- #
# Safety gate: require --yes-production to actually launch
# --------------------------------------------------------------------------- #
if [[ -z "$YES_PRODUCTION" ]]; then
    echo ""
    echo "========================================================"
    echo "  DRY PLAN COMPLETE (no shards launched)."
    echo "  Reference, validation, benchmark and shards are ready."
    echo "  Re-run with --yes-production to launch the $N_GPU shard(s)."
    echo "========================================================"
    exit 0
fi

# --------------------------------------------------------------------------- #
# Step 5 -- Launch shards (one per GPU), wait, collect exit codes
# --------------------------------------------------------------------------- #
echo "-- Step 5: launching $N_GPU shard(s)..."
PIDS=(); SHARDS=()
for ((i=0;i<N_GPU;i++)); do
    SHARD=$(printf "%s/shard_%03d.json" "$SHARD_DIR" "$i")
    LOG=$(printf "%s/shard_%03d.log" "$LOG_DIR" "$i")
    GPU=${GPU_LIST[$i]}
    echo "   shard_$(printf %03d $i) -> GPU $GPU  (log: $LOG)"
    CUDA_VISIBLE_DEVICES="$GPU" $PY scripts/run_gpu_shard.py \
        --config "$CONFIG" --stage "$STAGE" --shard "$SHARD" \
        > "$LOG" 2>&1 &
    PIDS+=($!); SHARDS+=("$SHARD")
done

FAIL=0
for ((i=0;i<${#PIDS[@]};i++)); do
    if ! wait "${PIDS[$i]}"; then
        echo "   !! shard ${SHARDS[$i]} exited non-zero (see logs)"; FAIL=1
    fi
done
echo "-- Step 5: all shards finished (fail flag=$FAIL)"

# --------------------------------------------------------------------------- #
# Step 6 -- Merge CSVs + config summaries
# --------------------------------------------------------------------------- #
echo "-- Step 6: merging shard CSVs + config summaries..."
$PY scripts/run_abf_fr_grid_torch.py --config "$CONFIG" --stage "$STAGE" --merge-only

# --------------------------------------------------------------------------- #
# Step 7 -- Report tables + plots
# --------------------------------------------------------------------------- #
echo "-- Step 7: report tables + figures..."
$PY scripts/make_report_tables.py --stage "$STAGE" --output-root "$OUTPUT_ROOT" || true
$PY scripts/plot_abf_fr_study.py --stage "$STAGE" --output-root "$OUTPUT_ROOT" || true

echo "========================================================"
echo "  GPU study complete (fail flag=$FAIL)."
echo "  finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Outputs : $STAGE_DIR/"
echo "    ${STAGE}_final_summary.csv, ${STAGE}_runs_long.csv, ..."
echo "    table_main_results.csv, best_configs.csv"
echo "    logs/  shards/  completed/  failed/"
echo "========================================================"
exit $FAIL
