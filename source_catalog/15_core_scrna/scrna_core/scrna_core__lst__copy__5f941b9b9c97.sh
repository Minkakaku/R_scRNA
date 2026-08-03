#!/bin/bash

OUTDIR="/data/outputdata/outputdata-sc/cellranger_out/FWSC20240705-14"
cd "$OUTDIR"
if [ -d "$OUTDIR" ]; then
    for dir in "$OUTDIR"/*; do
        if [ -d "$dir" ]; then
            # 定义子目录路径
            folder_name=$(basename "$dir")

            tmp="$OUTDIR/$folder_name/outs/"
            mkdir -p "$tmp"  # 确保目标目录存在

            # 复制 web_summary.html
            cp "$folder_name/outs/per_sample_outs/$folder_name/web_summary.html" "$tmp" 2>/dev/null || echo "$folder_name/outs/per_sample_outs/$folder_name/web_summary.html"

            # 复制 metrics_summary.csv
            cp "$folder_name/outs/per_sample_outs/$folder_name/metrics_summary.csv" "$tmp" 2>/dev/null || echo "未找到 $folder_name/outs/per_sample_outs/$folder_name/metrics_summary.csv"

            # 复制 filtered_feature_bc_matrix
            tmp1="$tmp/filtered_feature_bc_matrix/"
            mkdir -p "$tmp1"  # 确保目录存在
            cp "$folder_name/outs/per_sample_outs/$folder_name/count/sample_filtered_feature_bc_matrix/"* "$tmp1" 2>/dev/null || echo "$folder_name/outs/per_sample_outs/$folder_name/count/sample_filtered_feature_bc_matrix"
        fi
    done
fi