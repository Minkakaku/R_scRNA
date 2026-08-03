#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${SCRIPT_DIR}"
R_SCRIPT="${WORK_DIR}/DE_Analysis.r"
REPORT_SCRIPT="${WORK_DIR}/DE_Report.R"

CONFIG_FILE="${WORK_DIR}/compare_list.csv"
SPECIE=""
DATA_FILE=""
DATA_TYPE=""
GROUP_COL=""
CELLTYPE_COL=""
CPU=""
CONTAINER_PATH=""
REPORT_OUTDIR=""

usage() {
    cat <<EOF
Usage:
  bash DE_Analysis.sh [options]

Options:
  --config-file <file>      Compare list CSV. Optional, default: ${CONFIG_FILE}
  --data-file <file>        Input Seurat object file. Required.
  --data-type <type>        qs|rds|rdata|auto. Required.
  --species <type>          human|mouse. Required.
  --group-col <name>        Group column name. Required.
  --celltype-col <name>     Celltype column name. Required.
  --cpu <int>               CPU number. Required.
  --conda-path <dir>        Conda env prefix. Required.
  --output-dir <dir>        Output root for analysis/report. Required.
  -h, --help                Show this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config-file)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --data-file)
            DATA_FILE="$2"
            shift 2
            ;;
        --data-type)
            DATA_TYPE="$2"
            shift 2
            ;;
        --species)
            SPECIE="$2"
            shift 2
            ;;
        --group-col)
            GROUP_COL="$2"
            shift 2
            ;;
        --celltype-col)
            CELLTYPE_COL="$2"
            shift 2
            ;;
        --cpu)
            CPU="$2"
            shift 2
            ;;
        --conda-path)
            CONTAINER_PATH="$2"
            shift 2
            ;;
        --output-dir)
            REPORT_OUTDIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            exit 1
            ;;
    esac
done

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: $CONFIG_FILE not found!"
    exit 1
fi

if [ -z "$DATA_FILE" ]; then
    echo "Error: --data-file is required!"
    exit 1
fi

if [ ! -f "$DATA_FILE" ]; then
    echo "Error: $DATA_FILE not found!"
    exit 1
fi

if [ -z "$DATA_TYPE" ]; then
    echo "Error: --data-type is required!"
    exit 1
fi

if [ -z "$SPECIE" ]; then
    echo "Error: --species is required!"
    exit 1
fi

if [ -z "$GROUP_COL" ]; then
    echo "Error: --group-col is required!"
    exit 1
fi

if [ -z "$CELLTYPE_COL" ]; then
    echo "Error: --celltype-col is required!"
    exit 1
fi

if [ -z "$CPU" ]; then
    echo "Error: --cpu is required!"
    exit 1
fi

if [ -z "$CONTAINER_PATH" ]; then
    echo "Error: --conda-path is required!"
    exit 1
fi

if [ ! -d "$CONTAINER_PATH" ] && [ ! -f "$CONTAINER_PATH" ]; then
    echo "Error: Container path $CONTAINER_PATH not found!"
    exit 1
fi

if [ ! -f "$R_SCRIPT" ]; then
    echo "Error: Analysis script $R_SCRIPT not found!"
    exit 1
fi

if [ ! -f "$REPORT_SCRIPT" ]; then
    echo "Error: Report script $REPORT_SCRIPT not found!"
    exit 1
fi

if [ -z "$REPORT_OUTDIR" ]; then
    echo "Error: --output-dir is required!"
    exit 1
fi

mkdir -p "$REPORT_OUTDIR"
mkdir -p "${REPORT_OUTDIR}/log"

echo ">>> Start batch submission..."

tail -n +2 "$CONFIG_FILE" | tr -d '\r' | while IFS=, read -r ctrl_grp case_grp; do
    if [[ -z "$ctrl_grp" || -z "$case_grp" ]]; then
        continue
    fi

    echo ">>> Submitting: Case=$case_grp vs Control=$ctrl_grp"

    LOG_FILE="${REPORT_OUTDIR}/log/log_${case_grp}_vs_${ctrl_grp}.out"

    nohup conda run -p "$CONTAINER_PATH" \
        bash -lc "conda run -p '$CONTAINER_PATH' Rscript '$R_SCRIPT' '$case_grp' '$ctrl_grp' '$SPECIE' '$DATA_FILE' '$CPU' '$DATA_TYPE' '$GROUP_COL' '$CELLTYPE_COL' '$REPORT_OUTDIR' && conda run -p '$CONTAINER_PATH' Rscript '$REPORT_SCRIPT' '$case_grp' '$ctrl_grp' '$REPORT_OUTDIR'" \
        > "$LOG_FILE" 2>&1 &
done

echo ">>> All tasks have been submitted in background."
echo ">>> Use 'ps -ef | grep Rscript' to inspect progress."
