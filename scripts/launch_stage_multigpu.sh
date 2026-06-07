#!/usr/bin/env bash
# Launch one production stage sharded across GPUs 1-7 (GPU 0 is shared on this box).
# Each GPU runs one shard; runs are idempotent so re-launching resumes cleanly.
#
# Usage:
#   bash scripts/launch_stage_multigpu.sh <stage> [num_shards] [gpu_list_csv]
# Examples:
#   bash scripts/launch_stage_multigpu.sh pilot
#   bash scripts/launch_stage_multigpu.sh main 7 "1,2,3,4,5,6,7"
#   bash scripts/launch_stage_multigpu.sh difficulty_budget 5 "1,2,3,4,5"
#
# Logs: results/wca_production/logs/<stage>_shardK.log
# Waits for all shards, then prints a per-shard DONE summary.
set -u
cd /home/zheyuanlai/ABF-Fisher-Rao
source /home/zheyuanlai/miniconda3/etc/profile.d/conda.sh
conda activate abffr

STAGE="${1:?usage: launch_stage_multigpu.sh <stage> [num_shards] [gpu_csv]}"
GPUS_CSV="${3:-1,2,3,4,5,6,7}"
IFS=',' read -r -a GPUS <<< "$GPUS_CSV"
NUM_SHARDS="${2:-${#GPUS[@]}}"
CONFIG=configs/wca_production.yaml
LOGDIR=results/wca_production/logs
mkdir -p "$LOGDIR"

echo "[launch] stage=$STAGE num_shards=$NUM_SHARDS gpus=${GPUS_CSV}"
PIDS=()
for k in $(seq 0 $((NUM_SHARDS-1))); do
  gpu="${GPUS[$((k % ${#GPUS[@]}))]}"
  log="$LOGDIR/${STAGE}_shard${k}.log"
  echo "  shard $k -> GPU $gpu  ($log)"
  CUDA_VISIBLE_DEVICES="$gpu" nohup python -u scripts/run_wca_production.py \
    --config "$CONFIG" --stage "$STAGE" --shard "$k" --num-shards "$NUM_SHARDS" --verbose \
    > "$log" 2>&1 &
  PIDS+=("$!")
done
echo "[launch] launched ${#PIDS[@]} shards: ${PIDS[*]}"
echo "[launch] waiting..."
wait
echo "=== $STAGE DONE ==="
for k in $(seq 0 $((NUM_SHARDS-1))); do
  echo "### shard $k"; grep -E "\[run\] DONE|FAILED|NAN" "$LOGDIR/${STAGE}_shard${k}.log" | tail -3
done
