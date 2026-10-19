# 参考检查 -> hisat2 比对排序 -> (UMI 去重)


rule check_reference:
    """检查 hisat2 index 与 GTF 的染色体命名是否一致（chr1 vs 1），不一致时尽早报错。"""
    input:
        gtf="results/reference/genes.gtf",
    output:
        "results/qc/reference_check.txt",
    params:
        index=REF["hisat2_index"],
    log:
        "<logs>/check_reference.log",
    conda:
        ENV_MAIN
    shell:
        """
        hisat2-inspect -n {params.index} | awk '{{print $1}}' | sort -u > {output}.idx
        grep -v '^#' {input.gtf} | cut -f1 | sort -u > {output}.gtf
        n=$(comm -12 {output}.idx {output}.gtf | wc -l)
        {{
            echo "hisat2 index : {params.index}"
            echo "GTF          : {input.gtf}"
            echo "index 序列数 : $(wc -l < {output}.idx)"
            echo "GTF 序列数   : $(wc -l < {output}.gtf)"
            echo "共有序列数   : $n"
            echo "GTF 中不在 index 里的序列: $(comm -13 {output}.idx {output}.gtf | tr '\\n' ' ')"
        }} > {output}
        rm -f {output}.idx {output}.gtf
        cat {output} > {log}
        if [ "$n" -eq 0 ]; then
            echo "错误: hisat2 index 与 GTF 没有共同的染色体名（检查 chr 前缀是否一致）" | tee -a {log} >&2
            rm -f {output}
            exit 1
        fi
        """


def hisat2_reads(wc, input):
    if is_pe(wc.sample):
        return f"-1 {input.reads[0]} -2 {input.reads[1]}"
    return f"-U {input.reads[0]}"


rule hisat2:
    input:
        reads=trimmed_reads,
        check="results/qc/reference_check.txt",
    output:
        bam="results/03_align/{sample}.bam",
        bai="results/03_align/{sample}.bam.bai",
        summary="results/qc/hisat2/{sample}.summary.txt",  # --new-summary 格式，MultiQC 据此识别 hisat2
    params:
        index=REF["hisat2_index"],
        reads=hisat2_reads,
        strand=lambda wc: hisat2_strand(wc.sample),
        dta="--dta" if config["stringtie"]["enabled"] else "",
        extra=config["hisat2"]["params"],
        sort_mem=config["samtools_sort_mem"],
    threads: threads_for("align")
    log:
        hisat2="<logs>/hisat2/{sample}.log",
        sort="<logs>/hisat2/{sample}.sort.log",
    benchmark:
        "<benchmarks>/hisat2/{sample}.tsv"
    conda:
        ENV_MAIN
    shell:
        """
        hisat2 -p {threads} {params.extra} {params.dta} {params.strand} \
            --new-summary --summary-file {output.summary} \
            -x {params.index} {params.reads} 2> {log.hisat2} \
        | samtools sort -@ {threads} -m {params.sort_mem} -T {output.bam}.tmp -o {output.bam} - 2> {log.sort}
        samtools index -@ {threads} {output.bam} 2>> {log.sort}
        """


def umi_dedup_params(wc):
    u = config["umi"]
    args = "--extract-umi-method=read_id"
    if umi_in_comment(wc.sample):
        args += " --umi-separator='_'"  # umi_header 已改写为 name_UMI
    elif lib(wc.sample, "umi") == "header":
        args += f" --umi-separator='{u['separator']}'"
    if is_pe(wc.sample):
        args += f" --paired {u['dedup_extra_pe']}"
    return f"{args} {u['dedup_extra']}"


rule umi_dedup:
    input:
        bam="results/03_align/{sample}.bam",
        bai="results/03_align/{sample}.bam.bai",
    output:
        bam="results/03_align/{sample}.dedup.bam",
        bai="results/03_align/{sample}.dedup.bam.bai",
        stats="results/qc/umi/{sample}.dedup.log",  # umi_tools 统计（按比对记录计），MultiQC 使用
        counts="results/qc/umi/{sample}.read_counts.tsv",  # 去重前后的 reads 数（PE 为 pairs），summary 使用
    params:
        args=umi_dedup_params,
    threads: 1
    log:
        "<logs>/umi_dedup/{sample}.log",
    benchmark:
        "<benchmarks>/umi_dedup/{sample}.tsv"
    conda:
        ENV_MAIN
    shell:
        """
        umi_tools dedup {params.args} -I {input.bam} -S {output.bam} -L {output.stats} 2> {log}
        samtools index {output.bam} 2>> {log}
        # 多重比对的 read 在 BAM 中有多条记录，按 read 名去重计数；PE 两端同名，得到的是 pairs 数
        n() {{ samtools view -F 4 "$1" | cut -f1 | LC_ALL=C sort -u -S 1G | wc -l; }}
        printf "reads_before\treads_after\n%s\t%s\n" "$(n {input.bam})" "$(n {output.bam})" > {output.counts}
        """
