# 原始数据 -> (UMI extract / header 改写) -> trim_galore (-> 3' polyA 切除) ；raw / trimmed 各做一次 FastQC


rule link_raw:
    """把原始 fastq 统一链接为 {sample}_R1/R2.fastq.gz，后续文件名、报告中的样本名都与之一致。"""
    input:
        lambda wc: SAMPLES_DF.loc[wc.sample, "fq1" if wc.read == "R1" else "fq2"],
    output:
        "results/00_raw/{sample}_{read}.fastq.gz",
    localrule: True
    shell:
        "ln -sf {input:q} {output:q}"


# ------------------------------------------------------------
# UMI（umi=read 时）
# ------------------------------------------------------------
def umi_extract_params(wc):
    u = config["umi"]
    args = f"--extract-method={u['extract_method']} --bc-pattern='{u['pattern']}'"
    if is_pe(wc.sample) and u.get("pattern2"):
        args += f" --bc-pattern2='{u['pattern2']}'"
    return args


rule umi_extract_pe:
    wildcard_constraints:
        sample=regex_of([s for s in PE_SAMPLES if lib(s, "umi") == "read"]),
    input:
        r1="results/00_raw/{sample}_R1.fastq.gz",
        r2="results/00_raw/{sample}_R2.fastq.gz",
    output:
        r1=temp("results/01_umi/{sample}_R1.fastq.gz"),
        r2=temp("results/01_umi/{sample}_R2.fastq.gz"),
        stats="results/qc/umi/{sample}.extract.log",
    params:
        args=umi_extract_params,
    log:
        "<logs>/umi_extract/{sample}.log",
    conda:
        ENV_MAIN
    shell:
        """
        umi_tools extract {params.args} \
            -I {input.r1} --read2-in={input.r2} \
            -S {output.r1} --read2-out={output.r2} -L {output.stats} 2> {log}
        """


rule umi_extract_se:
    wildcard_constraints:
        sample=regex_of([s for s in SE_SAMPLES if lib(s, "umi") == "read"]),
    input:
        r1="results/00_raw/{sample}_R1.fastq.gz",
    output:
        r1=temp("results/01_umi/{sample}_R1.fastq.gz"),
        stats="results/qc/umi/{sample}.extract.log",
    params:
        args=umi_extract_params,
    log:
        "<logs>/umi_extract/{sample}.log",
    conda:
        ENV_MAIN
    shell:
        "umi_tools extract {params.args} -I {input.r1} -S {output.r1} -L {output.stats} 2> {log}"


# ------------------------------------------------------------
# UMI 在 header 注释中（umi=header 且 umi.header_field=comment）
# 比对软件会丢掉 read 名空格后的内容，先把 UMI 移到 read 名末尾：
#   @SRR123.1 V300...:2:CTCTCTAC:ACGCAGGCCC/1  ->  @SRR123.1_ACGCAGGCCC
# ------------------------------------------------------------
rule umi_header:
    wildcard_constraints:
        sample=regex_of([s for s in SAMPLES if umi_in_comment(s)]),
    input:
        "results/00_raw/{sample}_{read}.fastq.gz",
    output:
        temp("results/01_umi/{sample}_{read}.fastq.gz"),
    params:
        sep=config["umi"]["separator"],
    threads: threads_for("trim")
    log:
        "<logs>/umi_header/{sample}_{read}.log",
    conda:
        ENV_MAIN
    shell:
        """
        {{
        pigz -dc -p 2 {input} \
        | awk -v sep='{params.sep}' 'NR % 4 == 1 {{
            i = index($0, " ")
            if (i == 0) {{ print "第 " NR " 行 header 没有空格后的注释: " $0 > "/dev/stderr"; exit 1 }}
            c = substr($0, i + 1); sub(/\/[12]$/, "", c); n = split(c, f, sep)
            print substr($0, 1, i - 1) "_" f[n]; next
        }} {{ print }}' \
        | pigz -p {threads} > {output}
        }} 2> {log}
        """


# ------------------------------------------------------------
# trim_galore；protocol=3prime 时额外切除 5' 端并用 cutadapt 去 polyA
# ------------------------------------------------------------
def trim_args(wc):
    args = config["trim"]["params"]
    if lib(wc.sample, "protocol") == "3prime":
        p3 = config["trim"]["3prime"]
        if p3["clip_r1"]:
            args += f" --clip_R1 {p3['clip_r1']}"
        args += f" {p3['extra']}"
    return args


def polya_flag(wc):
    return int(lib(wc.sample, "protocol") == "3prime" and config["trim"]["3prime"]["polyA"])


rule trim_pe:
    wildcard_constraints:
        sample=regex_of(PE_SAMPLES),
    input:
        unpack(trim_input),
    output:
        r1=maybe_temp("results/02_trim/{sample}_R1.trimmed.fq.gz"),
        r2=maybe_temp("results/02_trim/{sample}_R2.trimmed.fq.gz"),
        report_r1="results/02_trim/reports/{sample}_R1.fastq.gz_trimming_report.txt",
        report_r2="results/02_trim/reports/{sample}_R2.fastq.gz_trimming_report.txt",
    params:
        args=trim_args,
        polya=polya_flag,
        minlen=config["trim"]["3prime"]["min_length"],
        tmp="results/02_trim/{sample}.tmp",
    threads: threads_for("trim")
    log:
        "<logs>/trim/{sample}.log",
    benchmark:
        "<benchmarks>/trim/{sample}.tsv"
    conda:
        ENV_MAIN
    shell:
        """
        rm -rf {params.tmp} && mkdir -p {params.tmp} results/02_trim/reports
        trim_galore {params.args} --cores {threads} --paired --basename {wildcards.sample} \
            -o {params.tmp} {input.r1} {input.r2} > {log} 2>&1
        if [ {params.polya} -eq 1 ]; then
            cutadapt -j {threads} -a 'A{{18}}' -m {params.minlen} \
                -o {output.r1} -p {output.r2} \
                {params.tmp}/{wildcards.sample}_val_1.fq.gz {params.tmp}/{wildcards.sample}_val_2.fq.gz \
                > results/02_trim/reports/{wildcards.sample}.polyA.cutadapt.log 2>> {log}
        else
            mv {params.tmp}/{wildcards.sample}_val_1.fq.gz {output.r1}
            mv {params.tmp}/{wildcards.sample}_val_2.fq.gz {output.r2}
        fi
        mv {params.tmp}/*_trimming_report.txt results/02_trim/reports/
        rm -rf {params.tmp}
        """


rule trim_se:
    wildcard_constraints:
        sample=regex_of(SE_SAMPLES),
    input:
        unpack(trim_input),
    output:
        r1=maybe_temp("results/02_trim/{sample}_R1.trimmed.fq.gz"),
        report_r1="results/02_trim/reports/{sample}_R1.fastq.gz_trimming_report.txt",
    params:
        args=trim_args,
        polya=polya_flag,
        minlen=config["trim"]["3prime"]["min_length"],
        tmp="results/02_trim/{sample}.tmp",
    threads: threads_for("trim")
    log:
        "<logs>/trim/{sample}.log",
    benchmark:
        "<benchmarks>/trim/{sample}.tsv"
    conda:
        ENV_MAIN
    shell:
        """
        rm -rf {params.tmp} && mkdir -p {params.tmp} results/02_trim/reports
        trim_galore {params.args} --cores {threads} --basename {wildcards.sample} \
            -o {params.tmp} {input.r1} > {log} 2>&1
        if [ {params.polya} -eq 1 ]; then
            cutadapt -j {threads} -a 'A{{18}}' -m {params.minlen} -o {output.r1} \
                {params.tmp}/{wildcards.sample}_trimmed.fq.gz \
                > results/02_trim/reports/{wildcards.sample}.polyA.cutadapt.log 2>> {log}
        else
            mv {params.tmp}/{wildcards.sample}_trimmed.fq.gz {output.r1}
        fi
        mv {params.tmp}/*_trimming_report.txt results/02_trim/reports/
        rm -rf {params.tmp}
        """


# ------------------------------------------------------------
# FastQC
# ------------------------------------------------------------
FASTQC_JAVA = f"--java {config['fastqc_java']}" if config.get("fastqc_java") else ""


rule fastqc_raw:
    input:
        "results/00_raw/{sample}_{read}.fastq.gz",
    output:
        html="results/qc/fastqc_raw/{sample}_{read}_fastqc.html",
        zip="results/qc/fastqc_raw/{sample}_{read}_fastqc.zip",
    params:
        java=FASTQC_JAVA,
    threads: threads_for("fastqc")
    log:
        "<logs>/fastqc_raw/{sample}_{read}.log",
    conda:
        ENV_MAIN
    shell:
        "fastqc {params.java} -t {threads} -o results/qc/fastqc_raw {input} > {log} 2>&1"


rule fastqc_trimmed:
    input:
        "results/02_trim/{sample}_{read}.trimmed.fq.gz",
    output:
        html="results/qc/fastqc_trimmed/{sample}_{read}.trimmed_fastqc.html",
        zip="results/qc/fastqc_trimmed/{sample}_{read}.trimmed_fastqc.zip",
    params:
        java=FASTQC_JAVA,
    threads: threads_for("fastqc")
    log:
        "<logs>/fastqc_trimmed/{sample}_{read}.log",
    conda:
        ENV_MAIN
    shell:
        "fastqc {params.java} -t {threads} -o results/qc/fastqc_trimmed {input} > {log} 2>&1"
