"""
每个样本一行的汇总表。reads 数按单条 read 计（PE 为 R1 + R2），所有 rate 都以 raw_reads 为分母（%）。
  raw_reads            trim 前（trim_galore 报告 Total reads processed，PE 两个报告相加）
  clean_reads          进入 hisat2 的 reads（trim / polyA 切除之后）
  mapped_reads         clean_reads x hisat2 overall alignment rate
  umi_dedup_reads      UMI 去重后的 reads（仅 umi != none；按 read 名计数，多重比对的 read 只算一次）
  clean_rate           clean_reads / raw_reads
  mapping_rate         mapped_reads / raw_reads
  umi_dedup_rate       umi_dedup_reads / raw_reads
  umi_duplicates_rate  (去重前比对上的 reads - umi_dedup_reads) / raw_reads
  assigned_rate        featureCounts Assigned / raw_reads（PE 的 fragment 数 x 2）
summary.tsv 的 reads 数带千位分隔符，MultiQC 用的表保持纯数字。
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")


def grab(pattern, text, cast=int, default=None):
    m = re.search(pattern, text)
    return cast(m.group(1).replace(",", "")) if m else default


def parse_trim(path):
    return grab(r"Total reads processed:\s+([\d,]+)", Path(path).read_text())


def parse_hisat2(path):
    """hisat2 --new-summary 格式；PE 的 Total 为 pairs，overall rate 按单条 read 计。"""
    t = Path(path).read_text()
    total = grab(r"Total (?:pairs|reads): (\d+)", t)
    rate = grab(r"Overall alignment rate: ([\d.]+)%", t, float)
    return total, rate


def pct(n, total):
    return 100 * n / total if n is not None and total else None


def parse_umi(path):
    r = pd.read_csv(path, sep="\t").iloc[0]
    return int(r["reads_before"]), int(r["reads_after"])


def parse_fc_summary(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    col = df.iloc[:, 0]
    return int(col.get("Assigned", 0))


rows = []
umi_counts = {Path(p).name.removesuffix(".read_counts.tsv"): p for p in snakemake.input.umi}
for i, sample in enumerate(snakemake.params.samples):
    meta = snakemake.params.meta[sample]
    div = 2 if meta["layout"] == "PE" else 1
    strand = pd.read_csv(snakemake.input.strand[i], sep="\t").iloc[0]
    raw = sum(parse_trim(fp) for fp in snakemake.params.trim[sample])
    clean, rate = parse_hisat2(snakemake.input.hisat2[i])
    clean = clean * div if clean is not None else None
    mapped = round(clean * rate / 100) if clean is not None and rate is not None else None
    assigned = parse_fc_summary(snakemake.input.counts[i] + ".summary") * div

    row = {
        "sample": sample,
        "layout": meta["layout"],
        "protocol": meta["protocol"],
        "umi": meta["umi"],
        "strand_provided": strand["provided"],
        "strand_inferred": strand["inferred"],
        "strand_used": strand["used"],
        "strand_status": strand["status"],
        "raw_reads": raw,
        "clean_reads": clean,
        "mapped_reads": mapped,
        "umi_dedup_reads": None,
        "clean_rate (%)": pct(clean, raw),
        "mapping_rate (%)": pct(mapped, raw),
        "umi_dedup_rate (%)": None,
        "umi_duplicates_rate (%)": None,
        "assigned_rate (%)": pct(assigned, raw),
    }
    if meta["umi"] != "none":
        n_in, n_out = (n * div for n in parse_umi(umi_counts[sample]))
        row["umi_dedup_reads"] = n_out
        row["umi_dedup_rate (%)"] = pct(n_out, raw)
        row["umi_duplicates_rate (%)"] = pct(n_in - n_out, raw)
    for k, v in row.items():
        if v is None and not k.startswith("umi_"):
            print(f"WARNING: {sample} 未能解析 {k}", file=sys.stderr)
    rows.append(row)

df = pd.DataFrame(rows)
if (df["umi"] == "none").all():
    df = df.drop(columns=[c for c in df.columns if c.startswith("umi_")])
count_cols = [c for c in df.columns if c.endswith("_reads")]
for c in df.columns:
    if c.endswith("(%)"):
        df[c] = df[c].astype(float).round(2)
    elif c in count_cols:
        df[c] = df[c].astype("Int64")

out = df.copy()
for c in count_cols:
    out[c] = out[c].map(lambda v: "" if pd.isna(v) else f"{int(v):,}")
out.to_csv(snakemake.output.tsv, sep="\t", index=False)

# MultiQC 自定义表格
with open(snakemake.output.mqc, "w") as f:
    f.write("# id: 'pipeline_summary'\n# section_name: 'Pipeline summary'\n"
            "# description: 'results/tables/summary.tsv（reads 按单条计，rate 以 raw_reads 为分母）'\n"
            "# plot_type: 'table'\n")
    df.to_csv(f, sep="\t", index=False)
