import gzip
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml


# ============================================================
# 环境与参考
# ============================================================
ENV_MAIN = config["envs"]["main"]
ENV_DEEPTOOLS = config["envs"]["deeptools"]

if config["genome"] not in config["references"]:
    raise ValueError(
        f"genome '{config['genome']}' 不在 references 中，可选: {list(config['references'])}"
    )
REF = config["references"][config["genome"]]

if config["repeats"]["enabled"]:
    _rep = config["repeats"]["gtf"]
    if not _rep or not Path(_rep).is_file():
        raise ValueError(
            f"repeats.enabled 为 true，需要提供 repeat GTF（例如由 RepeatMasker 注释转换的 GTF），"
            f"当前 repeats.gtf = '{_rep}'"
        )
    if Path(_rep).resolve() == Path(REF["gtf"]).resolve():
        print(
            f"警告: repeats.gtf 与基因注释是同一个文件（{_rep}），repeat 计数将与基因计数重复",
            file=sys.stderr,
        )

THREADS = config["threads"]


def threads_for(step):
    return THREADS.get(step, THREADS["default"])


# ============================================================
# 样本表
# ============================================================
LIB_OPTIONS = {
    "strandedness": ["auto", "unstranded", "forward", "reverse"],
    "protocol": ["full_length", "3prime"],
    "umi": ["none", "read", "header"],
}
# 链方向 -> featureCounts -s
STRAND_TO_S = {"unstranded": 0, "forward": 1, "reverse": 2}


def load_samples(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"样本表不存在: {path.resolve()}")
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#", index_col=False).fillna("")
    df.columns = df.columns.str.strip()
    for col in ("sample", "fq1"):
        if col not in df.columns:
            raise ValueError(f"样本表缺少必需列 '{col}': {path}")
    if "fq2" not in df.columns:
        df["fq2"] = ""
    df = df.apply(lambda c: c.str.strip())

    errors = []
    dup = df.loc[df["sample"].duplicated(), "sample"].tolist()
    if dup:
        errors.append(f"样本名重复: {dup}")
    for s in df["sample"]:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", s):
            errors.append(f"样本名只能包含字母、数字、. _ -（且不能以符号开头）: '{s}'")

    # 未填写的文库参数使用 config 默认值
    for col, choices in LIB_OPTIONS.items():
        default = config["library"][col]
        if col not in df.columns:
            df[col] = default
        df[col] = df[col].replace("", default)
        bad = df.loc[~df[col].isin(choices), ["sample", col]]
        for _, r in bad.iterrows():
            errors.append(f"{r['sample']}: {col}='{r[col]}' 无效，可选 {choices}")

    # fastq 路径：相对路径以样本表所在目录为基准
    for col in ("fq1", "fq2"):
        df[col] = [
            str((path.parent / p).resolve()) if p and not Path(p).is_absolute() else p
            for p in df[col]
        ]
    for _, r in df.iterrows():
        for col in ("fq1", "fq2"):
            if r[col] and not Path(r[col]).exists():
                errors.append(f"{r['sample']}: {col} 文件不存在: {r[col]}")
        if not r["fq1"]:
            errors.append(f"{r['sample']}: fq1 为空")

    if errors:
        raise ValueError("样本表检查未通过:\n  " + "\n  ".join(errors))
    return df.set_index("sample", drop=False)


SAMPLES_DF = load_samples(config["samples"])
SAMPLES = SAMPLES_DF["sample"].tolist()
PE_SAMPLES = [s for s in SAMPLES if SAMPLES_DF.loc[s, "fq2"]]
SE_SAMPLES = [s for s in SAMPLES if not SAMPLES_DF.loc[s, "fq2"]]


def umi_from_header(line, field, sep):
    """按 umi.header_field / umi.separator 从 fastq header 行中取出 UMI。"""
    name, _, comment = line[1:].rstrip("\n").partition(" ")
    text = comment if field == "comment" else name
    return re.sub(r"/[12]$", "", text).split(sep)[-1]


def check_umi_headers(n=1000):
    """抽查 umi=header 样本的前 n 条 reads，UMI 不合规时在启动阶段报错，避免静默去重失败。"""
    u = config["umi"]
    if u["header_field"] not in ("name", "comment"):
        raise ValueError(f"umi.header_field 只能是 name 或 comment，当前为 '{u['header_field']}'")
    errors = []
    for s in SAMPLES_DF.index[SAMPLES_DF["umi"] == "header"]:
        fq = SAMPLES_DF.loc[s, "fq1"]
        opener = gzip.open if fq.endswith(".gz") else open
        with opener(fq, "rt") as f:
            heads = [line for i, line in zip(range(4 * n), f) if i % 4 == 0]
        if u["header_field"] == "comment" and not all(" " in h for h in heads):
            errors.append(f"{s}: header 中没有空格后的注释，例: {heads[0].strip()}")
            continue
        umis = [umi_from_header(h, u["header_field"], u["separator"]) for h in heads]
        bad = [x for x in umis if not re.fullmatch(r"[ACGTN]+", x)]
        if bad or len({len(x) for x in umis}) > 1:
            errors.append(
                f"{s}: 取出的 UMI 不合规（应为等长的 ACGTN 序列），例: {heads[0].strip()} -> "
                f"'{umis[0]}'，长度 {sorted({len(x) for x in umis})}"
            )
    if errors:
        raise ValueError(
            "UMI 检查未通过，请检查 umi.header_field 与 umi.separator:\n  " + "\n  ".join(errors)
        )


check_umi_headers()


def regex_of(names):
    """把样本名列表转为 wildcard 约束；空列表返回永不匹配的正则。"""
    return "|".join(re.escape(n) for n in names) if names else "(?!x)x"


wildcard_constraints:
    sample=regex_of(SAMPLES),
    read="R1|R2",


def is_pe(sample):
    return sample in PE_SAMPLES


def lib(sample, key):
    return SAMPLES_DF.loc[sample, key]


def reads_of(sample):
    return ["R1", "R2"] if is_pe(sample) else ["R1"]


def samples_with(key, *values):
    return [s for s in SAMPLES if lib(s, key) in values]


# ============================================================
# 各步骤的输入
# ============================================================
def umi_in_comment(sample):
    return lib(sample, "umi") == "header" and config["umi"]["header_field"] == "comment"


def trim_input(wc):
    """UMI 在序列中时先 extract，UMI 在 header 注释中时先改写 read 名，否则直接用原始数据。"""
    if lib(wc.sample, "umi") == "read" or umi_in_comment(wc.sample):
        d = "results/01_umi"
    else:
        d = "results/00_raw"
    return {r.lower(): f"{d}/{wc.sample}_{r}.fastq.gz" for r in reads_of(wc.sample)}


def trimmed_reads(wc):
    return [f"results/02_trim/{wc.sample}_{r}.trimmed.fq.gz" for r in reads_of(wc.sample)]


def final_bam(wc):
    """定量用的 BAM：有 UMI 时用去重后的 BAM。"""
    if lib(wc.sample, "umi") == "none":
        return f"results/03_align/{wc.sample}.bam"
    return f"results/03_align/{wc.sample}.dedup.bam"


def hisat2_strand(sample):
    """auto 时比对阶段不加链参数（链方向在比对后推断）。"""
    strand = lib(sample, "strandedness")
    if strand in ("auto", "unstranded"):
        return ""
    if is_pe(sample):
        return "--rna-strandness " + ("FR" if strand == "forward" else "RF")
    return "--rna-strandness " + ("F" if strand == "forward" else "R")


def maybe_temp(path):
    return path if config["keep"]["trimmed_fastq"] else temp(path)


# ============================================================
# 目标文件
# ============================================================
def all_targets():
    t = []
    for s in SAMPLES:
        t += expand("results/qc/fastqc_raw/{s}_{r}_fastqc.zip", s=s, r=reads_of(s))
        t += expand("results/qc/fastqc_trimmed/{s}_{r}.trimmed_fastqc.zip", s=s, r=reads_of(s))
        if config["bigwig"]["enabled"]:
            t.append(f"results/07_bigwig/{s}.bw")
            if config["bigwig"]["split_strand"]:
                t += [f"results/07_bigwig/{s}.fwd.bw", f"results/07_bigwig/{s}.rev.bw"]
    t += [
        "results/tables/counts_exon.tsv",
        "results/tables/cpm_exon.tsv",
        "results/tables/counts_exon_intron.tsv",
        "results/tables/cpm_exon_intron.tsv",
        "results/tables/strandedness.tsv",
        "results/tables/summary.tsv",
        "results/qc/multiqc/multiqc_report.html",
    ]
    if samples_with("protocol", "full_length"):
        t += ["results/tables/tpm_exon.tsv", "results/tables/fpkm_exon.tsv"]
        if config["stringtie"]["enabled"]:
            t += ["results/tables/stringtie_tpm.tsv", "results/tables/stringtie_fpkm.tsv"]
    if config["repeats"]["enabled"]:
        t.append("results/tables/repeats_counts.tsv")
    return t


# ============================================================
# 运行记录：logs/run_<时间>/
# ============================================================
def pipeline_version(pipe):
    """流程当前的 commit；直接读 .git，运行时不需要 git。"""
    gd = pipe / ".git"
    try:
        head = (gd / "HEAD").read_text().strip()
        if not head.startswith("ref: "):
            return head
        ref = head[5:]
        if (gd / ref).exists():
            return f"{(gd / ref).read_text().strip()} ({ref})"
        for line in (gd / "packed-refs").read_text().splitlines():
            if line.endswith(" " + ref):
                return f"{line.split()[0]} ({ref})"
    except OSError:
        pass
    return "unknown"


def record_run(version_args):
    """onstart 调用：保存本次运行的配置、样本表、流程版本与软件版本。"""
    d = Path(RUN_DIR)
    d.mkdir(parents=True, exist_ok=True)
    with open(d / "config.yaml", "w") as f:
        yaml.safe_dump(dict(config), f, allow_unicode=True, sort_keys=False)
    shutil.copy(config["samples"], d / "samples.tsv")
    pipe = Path(workflow.basedir).parent
    (d / "pipeline.txt").write_text(
        f"pipeline: {pipe}\nversion:  {pipeline_version(pipe)}\ncommand:  {' '.join(sys.argv)}\n"
    )
    with open(d / "software_versions.tsv", "w") as f:
        subprocess.run(
            ["bash", str(pipe / "workflow/scripts/software_versions.sh"), *version_args],
            stdout=f, stderr=subprocess.DEVNULL,
            env={**os.environ, "SMK_PATH": sys.argv[0]},
        )


def finish_run(smk_log):
    """onsuccess/onerror 调用：合并本次各任务的 benchmark 为一张表，并保存 snakemake 主日志。"""
    d = Path(RUN_DIR)
    if not d.is_dir():  # dry-run
        return
    bdir = d / "benchmarks"
    tables = []
    for f in sorted(bdir.rglob("*.tsv")):
        df = pd.read_csv(f, sep="\t")
        df.insert(0, "job", f.stem)
        df.insert(0, "step", str(f.parent.relative_to(bdir)))
        tables.append(df)
    if tables:
        pd.concat(tables).to_csv(d / "benchmarks.tsv", sep="\t", index=False)
    shutil.rmtree(bdir, ignore_errors=True)
    if smk_log:
        shutil.copy(smk_log if isinstance(smk_log, str) else smk_log[0], d / "snakemake.log")
