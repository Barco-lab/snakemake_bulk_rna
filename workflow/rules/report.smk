# 汇总表、MultiQC、QC 报告


def trim_reports(sample):
    return [
        f"results/02_trim/reports/{sample}_{r}.fastq.gz_trimming_report.txt"
        for r in reads_of(sample)
    ]


def trim_report(s, r):
    return f"results/02_trim/reports/{s}_{r}.fastq.gz_trimming_report.txt"


rule summary:
    input:
        trim=[trim_report(s, r) for s in SAMPLES for r in reads_of(s)],
        hisat2=expand("results/qc/hisat2/{sample}.summary.txt", sample=SAMPLES),
        counts=expand("results/04_count/{sample}.featureCounts.txt", sample=SAMPLES),
        intron=expand("results/04_count/{sample}.intron.featureCounts.txt", sample=SAMPLES),
        strand=expand("results/04_count/strand/{sample}.tsv", sample=SAMPLES),
        umi=expand("results/qc/umi/{sample}.read_counts.tsv", sample=samples_with("umi", "read", "header")),
    output:
        "results/qc/summary.tsv",
    params:
        samples=SAMPLES,
        meta=SAMPLES_DF[["sample", "protocol", "umi"]]
        .assign(layout=["PE" if is_pe(s) else "SE" for s in SAMPLES])
        .to_dict("index"),
        trim={s: [trim_report(s, r) for r in reads_of(s)] for s in SAMPLES},
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
        summary="results/qc/summary.tsv",
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



def report_settings():
    """QC 报告「流程与参数」中的参数表：(中文名, 英文名, 中文取值, 英文取值)，取值中 `...` 显示为代码。"""
    c, ref = config, REF

    def count(col, zh=True):
        vc = SAMPLES_DF[col].value_counts().items()
        if zh:
            return "、".join(f"{v}（{n} 个样本）" for v, n in vc)
        return ", ".join(f"{v} ({n} samples)" for v, n in vc)

    u = c["umi"]
    zh, en = [], []
    if samples_with("umi", "header"):
        zh.append(f"header：UMI 在 fastq header 的 {u['header_field']} 字段，取按 `{u['separator']}` 分隔的最后一段")
        en.append(f"header: UMI is the last `{u['separator']}`-separated field of the fastq header {u['header_field']}")
    if samples_with("umi", "read"):
        pat = u["pattern"] + (f" / {u['pattern2']}" if u["pattern2"] else "")
        zh.append(f"read：`umi_tools extract --extract-method={u['extract_method']}`，pattern `{pat}`")
        en.append(f"read: `umi_tools extract --extract-method={u['extract_method']}`, pattern `{pat}`")
    if zh:
        extra = f"`{u['dedup_extra']}`" if u["dedup_extra"] else ""
        zh.append("比对后 `umi_tools dedup` 去重" + (f"，参数 {extra}" if extra else ""))
        en.append("deduplicated with `umi_tools dedup` after alignment" + (f", options {extra}" if extra else ""))
    t3 = c["trim"]["3prime"]
    exclude = c["gtf_filter"]["exclude_biotypes"]
    thr = c["strand_check"]["threshold"]
    clip = f"`--clip_R1 {t3['clip_r1']}`"
    return [
        ("样本数", "Samples", f"{len(SAMPLES)}（SE {len(SE_SAMPLES)}，PE {len(PE_SAMPLES)}）",
         f"{len(SAMPLES)} (SE {len(SE_SAMPLES)}, PE {len(PE_SAMPLES)})"),
        ("参考基因组", "Genome", c["genome"], c["genome"]),
        ("hisat2 index", "hisat2 index", f"`{ref['hisat2_index']}`", f"`{ref['hisat2_index']}`"),
        ("GTF", "GTF", f"`{ref['gtf']}`", f"`{ref['gtf']}`"),
        ("GTF 过滤", "GTF filter", "去掉 " + "、".join(exclude) if exclude else "不过滤",
         "remove " + ", ".join(exclude) if exclude else "none"),
        ("链方向（填写值）", "Strandedness (given)", count("strandedness") + f"；auto 时自动推断，阈值 {thr}",
         count("strandedness", False) + f"; inferred when auto, threshold {thr}"),
        ("文库类型", "Protocol", count("protocol") + "；full_length 输出 TPM/FPKM，3prime 只输出 counts/CPM",
         count("protocol", False) + "; TPM/FPKM for full_length only, counts/CPM for 3prime"),
        ("UMI", "UMI", count("umi"), count("umi", False)),
        ("UMI 处理", "UMI handling", "；".join(zh) or "不做 UMI 处理", "; ".join(en) or "none"),
        ("trim_galore", "trim_galore", f"`{c['trim']['params']}`", f"`{c['trim']['params']}`"),
        ("3' 端测序额外修剪", "Extra trimming for 3'-end", clip + ("，cutadapt 切除 polyA" if t3["polyA"] else ""),
         clip + (", polyA removed with cutadapt" if t3["polyA"] else "")),
        ("hisat2", "hisat2", f"`{c['hisat2']['params']}`", f"`{c['hisat2']['params']}`"),
        ("featureCounts", "featureCounts", f"`{c['featurecounts']['params']}`；另按 exon 优先计数内含子 reads",
         f"`{c['featurecounts']['params']}`; intronic reads counted exon-first"),
        ("bigwig", "bigwig", f"`{c['bigwig']['params']}`" if c["bigwig"]["enabled"] else "不输出",
         f"`{c['bigwig']['params']}`" if c["bigwig"]["enabled"] else "disabled"),
    ]


rule qc_report:
    """单文件 HTML 报告：概览、流程与参数、QC 汇总表（指标说明与参考值见 workflow/report/qc_metrics.yaml）。"""
    input:
        summary="results/qc/summary.tsv",
        strand="results/tables/strandedness.tsv",
        multiqc="results/qc/multiqc/multiqc_report.html",
        metrics=workflow.source_path("../report/qc_metrics.yaml"),
        flowchart=workflow.source_path("../../docs/images/bulk_rna_workflow.svg"),
    output:
        "results/qc/report.html",
    params:
        title=("bulk RNA-seq QC 报告", "bulk RNA-seq QC report"),
        pipeline="snakemake_bulk_rna",
        version=pipeline_version(Path(workflow.basedir).parent),
        settings=report_settings(),
        multiqc="multiqc/multiqc_report.html",  # 相对 report.html 的路径
    log:
        "<logs>/qc_report.log",
    localrule: True
    script:
        "../scripts/qc_report.py"
