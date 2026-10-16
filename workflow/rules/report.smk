# 汇总表、MultiQC


def trim_reports(sample):
    return [
        f"results/02_trim/reports/{sample}_{r}.fastq.gz_trimming_report.txt"
        for r in reads_of(sample)
    ]


rule summary:
    input:
        trim=[f"results/02_trim/reports/{s}_R1.fastq.gz_trimming_report.txt" for s in SAMPLES],
        hisat2=expand("results/qc/hisat2/{sample}.summary.txt", sample=SAMPLES),
        counts=expand("results/04_count/{sample}.featureCounts.txt", sample=SAMPLES),
        strand=expand("results/04_count/strand/{sample}.tsv", sample=SAMPLES),
        umi=expand("results/qc/umi/{sample}.read_counts.tsv", sample=samples_with("umi", "read", "header")),
    output:
        tsv="results/tables/summary.tsv",
        mqc="results/qc/summary_mqc.tsv",
    params:
        samples=SAMPLES,
        meta=SAMPLES_DF[["sample", "protocol", "umi"]]
        .assign(layout=["PE" if is_pe(s) else "SE" for s in SAMPLES])
        .to_dict("index"),
    log:
        "<logs>/summary.log",
    script:
        "../scripts/summary.py"


rule multiqc:
    input:
        fastqc=[
            f"results/qc/{d}/{s}_{r}{suffix}_fastqc.zip"
            for s in SAMPLES
            for r in reads_of(s)
            for d, suffix in (("fastqc_raw", ""), ("fastqc_trimmed", ".trimmed"))
        ],
        trim=[f for s in SAMPLES for f in trim_reports(s)],
        hisat2=expand("results/qc/hisat2/{sample}.summary.txt", sample=SAMPLES),
        counts=expand("results/04_count/{sample}.featureCounts.txt", sample=SAMPLES),
        umi=expand("results/qc/umi/{sample}.dedup.log", sample=samples_with("umi", "read", "header"))
        + expand("results/qc/umi/{sample}.extract.log", sample=samples_with("umi", "read")),
        summary="results/qc/summary_mqc.tsv",
        config=workflow.source_path("../assets/multiqc_config.yaml"),
    output:
        "results/qc/multiqc/multiqc_report.html",
    params:
        # featureCounts 的 .summary 不是规则输出，按路径传入
        extra=[f"results/04_count/{s}.featureCounts.txt.summary" for s in SAMPLES],
    log:
        "<logs>/multiqc.log",
    conda:
        ENV_MAIN
    shell:
        """
        multiqc -f -c {input.config} -o results/qc/multiqc -n multiqc_report.html \
            {input.fastqc} {input.trim} {input.hisat2} {input.umi} {input.summary} {params.extra} > {log} 2>&1
        """

