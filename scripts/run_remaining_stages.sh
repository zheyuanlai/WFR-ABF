#!/usr/bin/env bash
# Run the remaining stages sequentially across GPUs 1-7 (each shares the full pool,
# so they must run one stage at a time). Idempotent: re-running resumes.
# Usage: bash scripts/run_remaining_stages.sh
set -u
cd /home/zheyuanlai/ABF-Fisher-Rao
LOGDIR=results/wca_production/logs
GPUS="1,2,3,4,5,6,7"

for spec in "failure:7" "difficulty_budget:7" "difficulty_replicas:7" "difficulty_crowding:6"; do
  stage="${spec%%:*}"; nsh="${spec##*:}"
  echo "=== $(date +%H:%M:%S) launching $stage ($nsh shards) ==="
  bash scripts/launch_stage_multigpu.sh "$stage" "$nsh" "$GPUS" \
    > "$LOGDIR/${stage}_driver.log" 2>&1
  echo "=== $(date +%H:%M:%S) $stage finished ==="
  grep -hE "\[run\] DONE|FAILED" "$LOGDIR/${stage}_shard"*.log 2>/dev/null | sort | uniq -c
done
echo "=== ALL REMAINING STAGES DONE ==="
