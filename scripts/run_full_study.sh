#!/usr/bin/env bash
# run_full_study.sh — full production pipeline for the 2D ABF-FR study.
#
# Usage:
#   bash scripts/run_full_study.sh                          # production (tuning config)
#   bash scripts/run_full_study.sh --smoke                  # quick smoke check
#   bash scripts/run_full_study.sh --config path/to/cfg.yaml
#   bash scripts/run_full_study.sh --dry-run                # preview run list, no simulation
#   bash scripts/run_full_study.sh --skip-reference         # skip Step 1 (e.g. already done)
#   bash scripts/run_full_study.sh --tuning-only            # stop after tuning plots/tables
#
# Logs are written to results/two_dim_xi_x/run_full_study.log (tee'd to stdout).

set -euo pipefail

# --------------------------------------------------------------------------- #
# Defaults
# --------------------------------------------------------------------------- #
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

CONFIG="configs/two_dim_xi_x_tuning.yaml"
DRY_RUN=""
SKIP_REFERENCE=""
TUNING_ONLY=""
SMOKE=""

# --------------------------------------------------------------------------- #
# Argument parsing
# --------------------------------------------------------------------------- #
while [[ $# -gt 0 ]]; do
    case "$1" in
        --smoke)         SMOKE=1; CONFIG="configs/two_dim_xi_x_smoke.yaml"; shift ;;
        --config)        CONFIG="$2"; shift 2 ;;
        --dry-run)       DRY_RUN="--dry-run"; shift ;;
        --skip-reference) SKIP_REFERENCE=1; shift ;;
        --tuning-only)   TUNING_ONLY=1; shift ;;
        -h|--help)
            sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# --------------------------------------------------------------------------- #
# Logging setup
# --------------------------------------------------------------------------- #
LOG_DIR="results/two_dim_xi_x"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/run_full_study.log"

exec > >(tee -a "$LOG_FILE") 2>&1

STAGE_LABEL="production"
[[ -n "$SMOKE" ]] && STAGE_LABEL="smoke"

echo "========================================================"
echo "  ABF-FR 2D study — $STAGE_LABEL run"
echo "  config : $CONFIG"
echo "  started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "  log    : $LOG_FILE"
echo "========================================================"
echo ""

# --------------------------------------------------------------------------- #
# Helper
# --------------------------------------------------------------------------- #
run_step() {
    local label="$1"; shift
    echo "-------- $label — $(date '+%H:%M:%S') --------"
    local t0=$SECONDS
    "$@"
    local elapsed=$(( SECONDS - t0 ))
    echo "   done in ${elapsed}s"
    echo ""
}

# --------------------------------------------------------------------------- #
# Step 1 — Reference
# --------------------------------------------------------------------------- #
if [[ -z "$SKIP_REFERENCE" ]]; then
    run_step "Step 1: reference profiles and figures" \
        python3 scripts/run_reference_2d.py --config "$CONFIG"
else
    echo "-------- Step 1: reference — SKIPPED (--skip-reference) --------"
    echo ""
fi

# --------------------------------------------------------------------------- #
# Step 2 — Tuning grid
# --------------------------------------------------------------------------- #
STAGE="tuning"
[[ -n "$SMOKE" ]] && STAGE="smoke"

run_step "Step 2: ABF-FR grid (stage=$STAGE)" \
    python3 scripts/run_abf_fr_grid.py --config "$CONFIG" --stage "$STAGE" $DRY_RUN

# --------------------------------------------------------------------------- #
# Step 3 — Tuning plots and tables
# --------------------------------------------------------------------------- #
if [[ -z "$DRY_RUN" ]]; then
    run_step "Step 3a: tuning figures" \
        python3 scripts/plot_abf_fr_study.py --stage tuning

    run_step "Step 3b: tuning tables" \
        python3 scripts/make_report_tables.py --stage tuning
fi

if [[ -n "$TUNING_ONLY" ]]; then
    echo "========================================================"
    echo "  Stopped after tuning (--tuning-only)."
    echo "  Review results/two_dim_xi_x/tuning/ and figures_tuning/"
    echo "  then re-run without --tuning-only to continue to eval."
    echo "  finished: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "========================================================"
    exit 0
fi

# --------------------------------------------------------------------------- #
# Step 4 — Eval of selected configs
# --------------------------------------------------------------------------- #
if [[ -z "$DRY_RUN" ]]; then
    run_step "Step 4: eval of best configs (stage=eval)" \
        python3 scripts/run_abf_fr_grid.py --config "$CONFIG" --stage eval

    # --------------------------------------------------------------------------- #
    # Step 5 — Eval plots and tables
    # --------------------------------------------------------------------------- #
    run_step "Step 5a: eval figures" \
        python3 scripts/plot_abf_fr_study.py --stage eval

    run_step "Step 5b: eval tables" \
        python3 scripts/make_report_tables.py --stage eval
fi

echo "========================================================"
echo "  All steps complete."
echo "  finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "  Key outputs:"
echo "    tuning/best_configs.csv"
echo "    tuning/table_tuning_top_configs.csv"
echo "    eval/table_main_results.csv"
echo "    figures_tuning/  figures_eval/"
echo "========================================================"
