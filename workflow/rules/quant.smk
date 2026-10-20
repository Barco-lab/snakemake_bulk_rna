# GTF 过滤 -> 链方向检查 -> featureCounts（exon；内含子 exon 优先）-> 合并矩阵 (counts / CPM / TPM / FPKM)
# 可选: StringTie、repeat 计数


rule prepare_gtf:
    """按 gene_biotype 去掉不需要的基因（默认 miRNA/snoRNA/snRNA），支持 .gtf / .gtf.gz。"""
    input:
        REF["gtf"],
    output:
        "results/reference/genes.gtf",
    params:
        pattern=(
            'gene_(bio)?type "(' + "|".join(config["gtf_filter"]["exclude_biotypes"]) + ')"'
            if config["gtf_filter"]["exclude_biotypes"]
            else ""
        ),
    localrule: True
    shell:
        """
        if [[ {input} == *.gz ]]; then zcat {input}; else cat {input}; fi \
        | awk -v pat='{params.pattern}' 'pat == "" || $0 !~ pat' > {output}
        """


rule strand_check:
    """分别以 -s 1 / -s 2 计数，按 assigned 比例推断链方向；strandedness=auto 时采用推断值。"""
    input:
        bam=final_bam,
        gtf="results/reference/genes.gtf",
    output:
        tsv="results/04_count/strand/{sample}.tsv",
        s="results/04_count/strand/{sample}.s",
    params:
        provided=lambda wc: lib(wc.sample, "strandedness"),
        pe=lambda wc: is_pe(wc.sample),
        threshold=config["strand_check"]["threshold"],
        featurecounts=f"{ENV_DEEPTOOLS}/bin/featureCounts",
    threads: threads_for("count")
    log:
        "<logs>/strand_check/{sample}.log",
    script:
        "../scripts/strand_check.py"


rule featurecounts:
    input:
        bam=final_bam,
        gtf="results/reference/genes.gtf",
        s="results/04_count/strand/{sample}.s",
    output:
        "results/04_count/{sample}.featureCounts.txt",
    params:
        pe=lambda wc: "-p --countReadPairs" if is_pe(wc.sample) else "",
        args=config["featurecounts"]["params"],
        attrs=config["featurecounts"]["extra_attributes"],
    threads: threads_for("count")
    log:
        "<logs>/featurecounts/{sample}.log",
    benchmark:
        "<benchmarks>/featurecounts/{sample}.tsv"
    conda:
        ENV_DEEPTOOLS
    shell:
        """
        featureCounts -T {threads} -s $(cat {input.s}) {params.pe} {params.args} \
            --extraAttributes {params.attrs} \
            -a {input.gtf} -o {output} {input.bam} > {log} 2>&1
        """


def gene_body_args(args):
    """把 featurecounts.params 中的 -t 换成 gene（按基因全长计数）。"""
    args = re.sub(r"(^|\s)-t\s+\S+", " ", args)
    return f"{args} -t gene".strip()


rule featurecounts_intron:
    """
    内含子计数（参考 zUMIs，exon 优先）：先按 exon 计数并标记每条 read，
    未分配给 exon 的 reads 再按基因全长（-t gene）计数，一条 read 只算一次。
    """
    input:
        bam=final_bam,
        gtf="results/reference/genes.gtf",
        s="results/04_count/strand/{sample}.s",
    output:
        "results/04_count/{sample}.intron.featureCounts.txt",
    params:
        pe=lambda wc: "-p --countReadPairs" if is_pe(wc.sample) else "",
        args=config["featurecounts"]["params"],
        gene_args=gene_body_args(config["featurecounts"]["params"]),
        tmp="results/04_count/{sample}.intron_tmp",
    threads: threads_for("count")
    log:
        "<logs>/featurecounts_intron/{sample}.log",
    benchmark:
        "<benchmarks>/featurecounts_intron/{sample}.tsv"
    conda:
        ENV_DEEPTOOLS
    shell:
        """
        rm -rf {params.tmp} && mkdir -p {params.tmp}
        s=$(cat {input.s})
        featureCounts -T {threads} -s $s {params.pe} {params.args} -R BAM --Rpath {params.tmp} \
            -a {input.gtf} -o {params.tmp}/exon.txt {input.bam} > {log} 2>&1
        # featureCounts 只给分配成功的 read 加 XT；不能用 XS 判断（hisat2 自带 XS:A 链方向标签）
        samtools view -@ {threads} -b -e '!exists([XT])' -o {params.tmp}/not_exon.bam \
            {params.tmp}/*.featureCounts.bam 2>> {log}
        featureCounts -T {threads} -s $s {params.pe} {params.gene_args} \
            -a {input.gtf} -o {output} {params.tmp}/not_exon.bam >> {log} 2>&1
        rm -rf {params.tmp}
        """


rule merge_counts:
    input:
        counts=expand("results/04_count/{sample}.featureCounts.txt", sample=SAMPLES),
    output:
        counts="results/tables/counts_exon.tsv",
        cpm="results/tables/cpm_exon.tsv",
    params:
        samples=SAMPLES,
        mode="counts",
    log:
        "<logs>/merge/counts.log",
    script:
        "../scripts/merge_counts.py"


rule merge_counts_exon_intron:
    """exon + intron 计数（核 RNA 使用）；length 列为基因全长。"""
    input:
        counts=expand("results/04_count/{sample}.featureCounts.txt", sample=SAMPLES),
        add=expand("results/04_count/{sample}.intron.featureCounts.txt", sample=SAMPLES),
    output:
        counts="results/tables/counts_exon_intron.tsv",
        cpm="results/tables/cpm_exon_intron.tsv",
    params:
        samples=SAMPLES,
        mode="counts",
    log:
        "<logs>/merge/counts_exon_intron.log",
    script:
        "../scripts/merge_counts.py"


rule length_normalize:
    """TPM / FPKM 只对 full_length 样本计算（3' 端测序按长度归一化没有意义）。"""
    input:
        counts=expand(
            "results/04_count/{sample}.featureCounts.txt",
            sample=samples_with("protocol", "full_length"),
        ),
    output:
        tpm="results/tables/tpm_exon.tsv",
        fpkm="results/tables/fpkm_exon.tsv",
    params:
        samples=samples_with("protocol", "full_length"),
        mode="length",
    log:
        "<logs>/merge/tpm_fpkm.log",
    script:
        "../scripts/merge_counts.py"


rule merge_strandedness:
    input:
        expand("results/04_count/strand/{sample}.tsv", sample=SAMPLES),
    output:
        "results/tables/strandedness.tsv",
    localrule: True
    shell:
        "awk 'FNR > 1 || NR == 1' {input} > {output}"


# ------------------------------------------------------------
# StringTie（可选，仅 full_length）
# ------------------------------------------------------------
rule stringtie:
    input:
        bam=final_bam,
        gtf="results/reference/genes.gtf",
        s="results/04_count/strand/{sample}.s",
    output:
        gtf="results/05_stringtie/{sample}/{sample}.gtf",
        abund="results/05_stringtie/{sample}/gene_abund.tab",
    params:
        extra=config["stringtie"]["params"],
    threads: threads_for("count")
    log:
        "<logs>/stringtie/{sample}.log",
    benchmark:
        "<benchmarks>/stringtie/{sample}.tsv"
    conda:
        ENV_MAIN
    shell:
        """
        case $(cat {input.s}) in 1) strand=--fr ;; 2) strand=--rf ;; *) strand= ;; esac
        stringtie {input.bam} -p {threads} -e -G {input.gtf} $strand {params.extra} \
            -A {output.abund} -o {output.gtf} > {log} 2>&1
        """


rule merge_stringtie:
    input:
        expand(
            "results/05_stringtie/{sample}/gene_abund.tab",
            sample=samples_with("protocol", "full_length"),
        ),
    output:
        tpm="results/tables/stringtie_tpm.tsv",
        fpkm="results/tables/stringtie_fpkm.tsv",
    params:
        samples=samples_with("protocol", "full_length"),
    log:
        "<logs>/merge/stringtie.log",
    script:
        "../scripts/merge_stringtie.py"


# ------------------------------------------------------------
# Repeat 计数（可选，需要在 config 中提供 repeat GTF）
# ------------------------------------------------------------
rule featurecounts_repeats:
    input:
        bam=final_bam,
        s="results/04_count/strand/{sample}.s",
    output:
        "results/06_repeats/{sample}.featureCounts.txt",
    params:
        gtf=config["repeats"]["gtf"],
        pe=lambda wc: "-p --countReadPairs" if is_pe(wc.sample) else "",
        args=config["repeats"]["params"],
    threads: threads_for("count")
    log:
        "<logs>/featurecounts_repeats/{sample}.log",
    conda:
        ENV_DEEPTOOLS
    shell:
        """
        featureCounts -T {threads} -s $(cat {input.s}) {params.pe} {params.args} \
            -a {params.gtf} -o {output} {input.bam} > {log} 2>&1
        """


rule merge_repeats:
    input:
        counts=expand("results/06_repeats/{sample}.featureCounts.txt", sample=SAMPLES),
    output:
        counts="results/tables/repeats_counts.tsv",
    params:
        samples=SAMPLES,
        mode="counts_only",
    log:
        "<logs>/merge/repeats.log",
    script:
        "../scripts/merge_counts.py"
