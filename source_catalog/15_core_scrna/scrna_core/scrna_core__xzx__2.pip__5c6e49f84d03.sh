#!/bin/bash
# 2.pip.sh

# 获取参数
meta_file="$1"
method="$2"
group="$3"
outdir="$4"
method_color="$5"
group_color="$6"

# 获取当前用户信息
USER_ID=$(id -u)
GROUP_ID=$(id -g)
USER_NAME=$(id -un)

# 获取绝对路径
OUTPUT_DIR=$(realpath "$outdir")
PROP_DIR=$(realpath "/mnt/NexData_file/share/scRNA-tools/Proption")

# 调试信息
echo "=== Debug Info ==="
echo "User: $USER_NAME ($USER_ID:$GROUP_ID)"
echo "Meta file: $meta_file"
echo "Method: $method"
echo "Group: $group"
echo "Output dir: $outdir ($OUTPUT_DIR)"
echo "Method color: $method_color"
echo "Group color: $group_color"
echo "=================="

docker run --name proportion_wzb --rm \
  -v "$OUTPUT_DIR:$OUTPUT_DIR" \
  -v "$PROP_DIR:$PROP_DIR" \
#  -u $USER_NAME:$USER_NAME \
  wangzhenbo/monocle:20240909 /bin/bash -c "
    echo '=== Inside Container ==='
    echo 'Checking file existence:'
    ls -l \"$meta_file\" \"$method_color\" \"$group_color\"
    echo '======================='
    
    # 设置文件权限
    echo $(whoami)
    
    # 执行R脚本
    
    Rscript \"$PROP_DIR/Rscript/cell_proportion_v1.R\" \
        -m \"$meta_file\" \
        -c \"$method\" \
        -g \"$group\" \
        -o \"$outdir\" \
        -a \"$method_color\" \
        -b \"$group_color\"

  "












# #!/bin/bash
# # run_proportion.sh

# # 获取参数
# meta_file=$1
# method=$2
# group=$3
# outdir=$4
# method_color=$5
# group_color=$6

# # 获取当前用户信息
# USER_ID=$(id -u)
# GROUP_ID=$(id -g)
# USER_NAME=$(id -un)  # 修正了这里的语法错误

# # 获取绝对路径
# OUTPUT_DIR=$(realpath "$outdir")
# PROP_DIR=$(realpath "/mnt/NexData_file/share/scRNA-tools/Proption")

# docker run --name proportion_wzb --rm \
#   -v "$OUTPUT_DIR:$OUTPUT_DIR" \
#   -v "$PROP_DIR:$PROP_DIR" \
#   -v /home/wzb:/home/wzb \
#   -e USER_ID=$USER_ID \
#   -e GROUP_ID=$GROUP_ID \
#   wangzhenbo/monocle:20240909 /bin/bash -c "
#     # 确保用户存在
#     if ! id -u $USER_NAME >/dev/null 2>&1; then
#       groupadd -g $GROUP_ID $USER_NAME || true
#       useradd -u $USER_ID -g $GROUP_ID -m $USER_NAME || true
#     fi
    
#     # 以主机用户身份运行命令
#     su - $USER_NAME -c \"
#       Rscript $PROP_DIR/Rscript/cell_proportion_v1.R \
#         -m '$meta_file' \
#         -c '$method' \
#         -g '$group' \
#         -o '$outdir' \
#         -a '$method_color' \
#         -b '$group_color'
#     \"
#   "