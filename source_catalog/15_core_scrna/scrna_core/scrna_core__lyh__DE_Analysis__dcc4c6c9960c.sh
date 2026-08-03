#!/bin/bash

# ================= 配置区域 =================
# R 脚本路径 (建议使用绝对路径，防止容器内路径混淆)
# 假设当前目录在 /public23/... 下，最好获取当前脚本的绝对路径
WORK_DIR=$(pwd)
R_SCRIPT="${WORK_DIR}/DE_Analysis.r"

# 配置文件路径
CONFIG_FILE="compare_list.csv"

# 参数设置
SPECIE="human"
QSDATA="/home/lyh/sc_RNA_seq/Jing/cellanno_T/anno.qs"
CPU=8

# 容器配置
CONTAINER_PATH="/mnt/NexData_file/home/lyh/Software/sif/my_sandbox"
BIND_PATHS="/mnt,/home"
# ===========================================

# 检查必要文件
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: $CONFIG_FILE not found!"
    exit 1
fi

if [ ! -d "$CONTAINER_PATH" ] && [ ! -f "$CONTAINER_PATH" ]; then
    echo "Error: Container path $CONTAINER_PATH not found!"
    exit 1
fi

echo ">>> 开始批量投递任务 (Container Mode)..."

# 读取 CSV (跳过第一行 header)
tail -n +2 "$CONFIG_FILE" | tr -d '\r' | while IFS=, read -r ctrl_grp case_grp; do
    
    # 跳过空行
    if [[ -z "$ctrl_grp" || -z "$case_grp" ]]; then continue; fi

    echo ">>> Submitting: Case=$case_grp vs Control=$ctrl_grp"

    # 定义日志文件
    LOG_FILE="${WORK_DIR}/log_${case_grp}_vs_${ctrl_grp}.out"

    # === 核心修改：使用容器执行 ===
    # singularity exec -B <挂载路径> <容器路径> <命令>
    # nohup 放在最外层以保持后台运行
    nohup singularity exec \
        -B "$BIND_PATHS" \
        "$CONTAINER_PATH" \
        Rscript "$R_SCRIPT" "$case_grp" "$ctrl_grp" "$SPECIE" "$QSDATA" "$CPU" \
        > "$LOG_FILE" 2>&1 &

done

echo ">>> 所有任务已投递到后台 (容器环境)。"
echo ">>> 使用 'ps -ef | grep Rscript' 或 'ps -ef | grep singularity' 查看进程。"