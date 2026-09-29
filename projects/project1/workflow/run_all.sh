#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$project_dir/logs"
for step in 00_setup 01_qc; do
  Rscript "$project_dir/workflow/$step.R" 2>&1 | tee "$project_dir/logs/$step.log"
done
printf '%s\n' '检查 config/qc_thresholds.csv 和 samples.tsv 后再运行过滤与聚类；注释前再次人工检查。'
