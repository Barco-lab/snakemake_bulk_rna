"""
每个样本一行的汇总表。reads 数按单条 read 计（PE 为 R1 + R2）。
rate 以进入该步骤的 reads 为分母（与各软件报告、MultiQC 及公开参考值一致），分母写在表头中：
  clean_rate 以 raw，mapping_rate 以 clean，UMI 去重以去重前比对上的 reads，
  assigned / intron 以进入 featureCounts 的 reads（UMI 样本为去重后，其余为 mapped）为分母。
输出纯数字，MultiQC 与 qc_report 读取；各列含义、参考值见 workflow/report/qc_metrics.yaml。
"""
import re
import sys
from pathlib import Path

import pandas as pd


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
    intron = parse_fc_summary(snakemake.input.intron[i] + ".summary") * div

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
        "clean_rate (% of raw)": pct(clean, raw),
        "mapping_rate (% of clean)": pct(mapped, clean),
        "umi_dedup_rate (% of mapped)": None,
        "umi_duplicates_rate (% of mapped)": None,
    }
    counted = mapped  # 进入 featureCounts 的 reads
    if meta["umi"] != "none":
        n_in, n_out = (n * div for n in parse_umi(umi_counts[sample]))
        row["umi_dedup_reads"] = n_out
        row["umi_dedup_rate (% of mapped)"] = pct(n_out, n_in)
        row["umi_duplicates_rate (% of mapped)"] = pct(n_in - n_out, n_in)
        counted = n_out
    row["assigned_rate (% of counted)"] = pct(assigned, counted)
    row["intron_rate (% of counted)"] = pct(intron, counted)
    for k, v in row.items():
        if v is None and not k.startswith("umi_"):
            print(f"WARNING: {sample} 未能解析 {k}", file=sys.stderr)
    rows.append(row)

df = pd.DataFrame(rows)
if (df["umi"] == "none").all():
    df = df.drop(columns=[c for c in df.columns if c.startswith("umi_")])
count_cols = [c for c in df.columns if c.endswith("_reads")]
for c in df.columns:
    if "(%" in c:
        df[c] = df[c].astype(float).round(2)
    elif c in count_cols:
        df[c] = df[c].astype("Int64")

df.to_csv(snakemake.output[0], sep="\t", index=False)
