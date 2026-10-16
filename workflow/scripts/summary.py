"""
每个样本一行的汇总表。reads 数按"片段"计：PE 为 read pairs，SE 为 reads，比对率分子分母口径一致。
  raw_reads      trim 前（trim_galore 报告 Total reads processed）
  clean_reads    进入 hisat2 的 reads（trim / polyA 切除之后）
  mapping_rate   hisat2 overall alignment rate
  unique_rate    唯一比对（PE 为 concordantly exactly 1 time）
  umi_dedup_*    UMI 去重（仅 umi != none）：按 read 名计数（PE 为 pairs），多重比对的 read 只算一次
                 umi_dup_rate = 1 - 去重后 reads / 去重前比对上的 reads
  assigned_*     featureCounts Assigned
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
    """hisat2 --new-summary 格式。"""
    t = Path(path).read_text()
    total = grab(r"Total (?:pairs|reads): (\d+)", t)
    rate = grab(r"Overall alignment rate: ([\d.]+)%", t, float)
    uniq = grab(r"Aligned concordantly 1 time: (\d+)", t)
    if uniq is None:
        uniq = grab(r"Aligned 1 time: (\d+)", t)
    return total, rate, uniq


def parse_umi(path):
    r = pd.read_csv(path, sep="\t").iloc[0]
    return int(r["reads_before"]), int(r["reads_after"])


def parse_fc_summary(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    col = df.iloc[:, 0]
    return int(col.get("Assigned", 0)), int(col.sum())


rows = []
umi_counts = {Path(p).name.removesuffix(".read_counts.tsv"): p for p in snakemake.input.umi}
for i, sample in enumerate(snakemake.params.samples):
    meta = snakemake.params.meta[sample]
    strand = pd.read_csv(snakemake.input.strand[i], sep="\t").iloc[0]
    raw = parse_trim(snakemake.input.trim[i])
    clean, rate, uniq = parse_hisat2(snakemake.input.hisat2[i])
    assigned, fc_total = parse_fc_summary(snakemake.input.counts[i] + ".summary")

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
        "clean_ratio": clean / raw if raw and clean is not None else None,
        "mapping_rate": rate / 100 if rate is not None else None,
        "unique_rate": uniq / clean if clean and uniq is not None else None,
        "umi_dedup_reads": None,
        "umi_dup_rate": None,
        "assigned": assigned,
        "assigned_rate": assigned / fc_total if fc_total else None,
    }
    if meta["umi"] != "none":
        n_in, n_out = parse_umi(umi_counts[sample])
        row["umi_dedup_reads"] = n_out
        row["umi_dup_rate"] = 1 - n_out / n_in if n_in else None
    for k, v in row.items():
        if v is None and not k.startswith("umi_"):
            print(f"WARNING: {sample} 未能解析 {k}", file=sys.stderr)
    rows.append(row)

df = pd.DataFrame(rows)
if (df["umi"] == "none").all():
    df = df.drop(columns=["umi_dedup_reads", "umi_dup_rate"])
for c in df.columns:
    if c.endswith(("_rate", "_ratio")):
        df[c] = df[c].astype(float).round(4)
    elif c.endswith(("_reads", "assigned")):
        df[c] = df[c].astype("Int64")
df.to_csv(snakemake.output.tsv, sep="\t", index=False)

# MultiQC 自定义表格
with open(snakemake.output.mqc, "w") as f:
    f.write("# id: 'pipeline_summary'\n# section_name: 'Pipeline summary'\n"
            "# description: 'results/tables/summary.tsv（PE 按 read pairs 计）'\n"
            "# plot_type: 'table'\n")
    df.to_csv(f, sep="\t", index=False)
