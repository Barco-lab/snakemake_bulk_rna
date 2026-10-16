#!/usr/bin/env bash
# 输出流程用到的软件版本与实际路径（流程启动时由 onstart 调用）
#   software_versions.sh <main 环境> <omics 环境（deeptools、featureCounts）> [fastqc 的 --java 参数]
#   环境变量 SMK_PATH：正在运行的 snakemake
main=$1 dt=$2 java=${3:-}
export PATH=$main/bin:$PATH

v() { printf "%s\t%s\t%s\n" "$1" "$2" "$(command -v "$3" || echo "$3")"; }
printf "tool\tversion\tpath\n"
v snakemake     "$("$SMK_PATH" --version)" "$SMK_PATH"
v trim_galore   "$(trim_galore --version | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)" trim_galore
v cutadapt      "$(cutadapt --version)" cutadapt
v fastqc        "$(fastqc $java --version 2>/dev/null | grep -oE 'v[0-9.]+')" fastqc
v umi_tools     "$(umi_tools --version | grep -oE '[0-9.]+$')" umi_tools
v hisat2        "$(hisat2 --version | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')" hisat2
v samtools      "$(samtools --version | head -1 | cut -d' ' -f2)" samtools
v featureCounts "$("$dt/bin/featureCounts" -v 2>&1 | grep -oE 'v[0-9.]+')" "$dt/bin/featureCounts"
v stringtie     "$(stringtie --version)" stringtie
v multiqc       "$(multiqc --version | grep -oE '[0-9.]+$')" multiqc
v deeptools     "$("$dt/bin/bamCoverage" --version | grep -oE '[0-9.]+$')" "$dt/bin/bamCoverage"
