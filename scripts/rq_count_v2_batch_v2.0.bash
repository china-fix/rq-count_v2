#!/bin/bash

# 显示帮助文档
usage() {
    cat <<EOF
用法: $0 [选项]
对一个文件夹中的多个a样式文件，使用同一个b文件进行处理，并将结果输出到新的文件夹中。

选项:
  --input_dir DIR      输入文件夹，包含多个a样式的文件
  --index_file FILE    输入文件b，包含要定位的关键字段
  --output_dir DIR     输出文件夹，用于存放生成的c文件
  -t, --threads N      并行处理的线程数（默认为所有可用CPU核心）
  --plus VALUE         对b文件的第二列进行递增预处理（默认0，表示不进行处理）
  -h, --help           显示此帮助信息

示例:

a样式文件示例 (tab分隔):
  CM001153.1    1    0
  CM001153.1    2    0

b文件示例 (tab分隔):
  name1    1
  name2    3

输出文件c示例:
  name1    1    CM001153.1    1    0
  name2    3    CM001153.1    3    0
EOF
}

# 默认参数
threads=$(nproc) # 默认线程数为CPU核心数
plus=0           # 默认不进行递增处理

# 解析命令行参数
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_dir) input_dir="$2"; shift ;;
        --index_file) index_file="$2"; shift ;;
        --output_dir) output_dir="$2"; shift ;;
        -t|--threads) threads="$2"; shift ;;
        --plus) plus="$2"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "未知选项: $1"; usage; exit 1 ;;
    esac
    shift
done

# 检查必需的参数是否提供
if [ -z "$input_dir" ] || [ -z "$index_file" ] || [ -z "$output_dir" ]; then
    echo "错误: 必须提供 --input_dir, --index_file, 和 --output_dir 参数。"
    usage
    exit 1
fi

# 检查输入文件夹和文件是否存在
if [ ! -d "$input_dir" ]; then
    echo "错误: 输入文件夹 $input_dir 不存在！"
    exit 1
fi

if [ ! -f "$index_file" ]; then
    echo "错误: 输入文件 $index_file 不存在！"
    exit 1
fi

# 创建输出文件夹（如果不存在）
mkdir -p "$output_dir"

# 如果 --plus 参数大于 0，进行 b 文件的预处理
if [ "$plus" -gt 0 ]; then
    index_file_p="${index_file}_p${plus}"
    > "$index_file_p" # 清空或创建文件

    # 原始 b 文件内容直接复制到新文件
    cat "$index_file" >> "$index_file_p"

    # 对 b 文件第二列进行递增处理
    for ((i=1; i<=plus; i++)); do
        awk -v increment="$i" 'BEGIN {OFS="\t"} {print $1, $2 + increment}' "$index_file" >> "$index_file_p"
    done

    echo "已生成预处理文件: $index_file_p"
    index_file="$index_file_p" # 更新为新生成的文件
fi

# 定义单个文件的处理函数
process_file() {
    file_a="$1"
    index_file="$2"
    output_dir="$3"
    filename=$(basename "$file_a")
    output_file="$output_dir/$filename"

    # 创建输出文件
    > "$output_file"

    # 处理文件 b 的每一行
    while IFS=$'\t' read -r name number; do
        # 查找文件 a 中对应的行
        match=$(awk -v num="$number" '$2 == num' "$file_a")
        if [ -n "$match" ]; then
            # 如果找到对应的行，将其内容添加到文件 b 的行后面
            echo -e "$name\t$number\t$match" >> "$output_file"
        else
            # 如果没有找到，保持原样
            echo -e "$name\t$number" >> "$output_file"
        fi
    done < "$index_file"

    echo "已处理文件: $file_a -> $output_file"
}

export -f process_file

# 使用 GNU Parallel 并行处理
find "$input_dir" -type f | parallel -j "$threads" process_file {} "$index_file" "$output_dir"

echo "所有文件处理完成，结果已保存到 $output_dir"
