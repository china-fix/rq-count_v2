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
  -h, --help           显示此帮助信息

示例:
  ./process_multiple_files_parallel.sh --input_dir input_dir --index_file b.txt --output_dir output_dir --threads 4
EOF
}

# 默认线程数为系统核心数
threads=$(nproc)

# 解析命令行参数
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_dir) input_dir="$2"; shift ;;
        --index_file) index_file="$2"; shift ;;
        --output_dir) output_dir="$2"; shift ;;
        -t|--threads) threads="$2"; shift ;;
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
