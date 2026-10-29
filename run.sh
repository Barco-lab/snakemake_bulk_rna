#!/usr/bin/env bash
# bulk RNA-seq 流程启动脚本（本地运行）
#
#   bash run.sh -c <config.yaml> -i <fq_dir> [-o <outdir>] [-j 32] [-n] [-y] [-- <其它 snakemake 参数>]
#
# -c 指定的配置文件不存在时，会从默认配置生成一份模板，修改后再次运行即可。
# 首次运行会在 <outdir> 中生成 samples.tsv，确认无误后再次运行即可。
set -euo pipefail

PIPE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SMK_ENV=${SMK_ENV:-$HOME/.conda/envs/omics}   # 提供 snakemake（>= 9.12）的环境

usage() {
    cat <<EOF
用法: bash run.sh -c <config.yaml> [选项] [-- <snakemake 参数>]

  -c FILE  项目配置文件（必填）。文件不存在时从默认配置生成模板后退出
           可以只写与默认配置（$PIPE_DIR/config/config.yaml）不同的参数
  -i DIR   fastq 目录（.fastq.gz / .fq.gz；_1/_2、_R1/_R2 自动识别为双端）
  -o DIR   输出目录（结果、日志、samples.tsv 都在这里；默认 ./rnaseq_out）
  -s FILE  使用已有的样本表（不扫描 -i）
  -j INT   最多同时使用的 CPU 核数（默认 32）
  -n       只预览（dry-run），不实际运行
  -y       跳过运行前的确认
  -h       显示帮助

示例:
  bash run.sh -c project1/config.yaml                             # 生成配置模板
  bash run.sh -c project1/config.yaml -i rawdata -o project1 -n   # 生成样本表并预览
  nohup bash run.sh -c project1/config.yaml -o project1 -j 48 -y > project1/run.log 2>&1 &
  bash run.sh -c project1/config.yaml -o project1 -n -- --until featurecounts   # 只跑到某一步
EOF
}

die() { echo "错误: $*" >&2; exit 1; }

config="" fq_dir="" outdir="" samples="" cores=32 dry=0 yes=0
while getopts "c:i:o:s:j:nyh" opt; do
    case $opt in
        c) config=$OPTARG ;;
        i) fq_dir=$OPTARG ;;
        o) outdir=$OPTARG ;;
        s) samples=$OPTARG ;;
        j) cores=$OPTARG ;;
        n) dry=1 ;;
        y) yes=1 ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done
shift $((OPTIND - 1))
[ "${1:-}" = "--" ] && shift
extra=("$@")

# ---- config ----
[ -n "$config" ] || { usage; die "必须用 -c 指定配置文件"; }
if [ ! -f "$config" ]; then
    mkdir -p "$(dirname "$config")"
    cp "$PIPE_DIR/config/config.yaml" "$config"
    echo ">> 已生成配置模板 $config"
    echo "   请修改其中的物种、文库类型等参数（不需要改的参数可以删掉，会使用默认值），然后再次运行。"
    exit 0
fi
config=$(realpath "$config")

[ -x "$SMK_ENV/bin/snakemake" ] || die "找不到 $SMK_ENV/bin/snakemake（可通过 SMK_ENV=<环境路径> 指定）"
# 注意：不要把 $SMK_ENV/bin 加到 PATH 前面，否则会挡住各规则 conda 激活的环境

outdir=${outdir:-rnaseq_out}
[ -d "$outdir" ] || echo ">> 创建输出目录 $outdir"
mkdir -p "$outdir"
outdir=$(realpath "$outdir")

# ---- samples.tsv ----
if [ -n "$samples" ]; then
    [ -f "$samples" ] || die "样本表不存在: $samples"
    [ "$(realpath "$samples")" = "$outdir/samples.tsv" ] || cp "$samples" "$outdir/samples.tsv"
elif [ -f "$outdir/samples.tsv" ]; then
    [ -n "$fq_dir" ] && echo ">> $outdir/samples.tsv 已存在，沿用该文件（如需重新扫描请先删除它）"
else
    [ -n "$fq_dir" ] || die "首次运行需要 -i <fq_dir> 或 -s <samples.tsv>"
    [ -d "$fq_dir" ] || die "目录不存在: $fq_dir"
    "$SMK_ENV/bin/python" "$PIPE_DIR/workflow/scripts/make_samplesheet.py" -i "$fq_dir" -o "$outdir/samples.tsv" \
        --extra-cols strandedness,protocol,umi
fi
echo ">> 配置 $config"
echo ">> 样本表 $outdir/samples.tsv:"
column -t -s $'\t' "$outdir/samples.tsv" | sed 's/^/   /'

# 先读默认配置，再用 -c 的配置覆盖（嵌套的键逐层合并）
smk=("$SMK_ENV/bin/snakemake"
    --snakefile "$PIPE_DIR/workflow/Snakefile"
    --directory "$outdir"
    --configfile "$config"
    --workflow-profile "$PIPE_DIR/workflow/profiles/default"
    --cores "$cores")

echo ">> 预览将要运行的任务..."
# 额外参数放在 --quiet rules 之前，否则目标文件会被当作 --quiet 的取值
"${smk[@]}" ${extra[@]+"${extra[@]}"} --dry-run --quiet rules
[ "$dry" -eq 1 ] && exit 0

if [ "$yes" -ne 1 ]; then
    read -r -p "确认开始运行？(y/n): " ans
    [[ "$ans" =~ ^[Yy]$ ]] || { echo "已取消"; exit 0; }
fi

echo ">> 开始运行 $(date '+%F %T')"
"${smk[@]}" ${extra[@]+"${extra[@]}"}
echo ">> 结束 $(date '+%F %T')，结果在 $outdir/results"
