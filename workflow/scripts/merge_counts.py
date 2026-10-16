"""
合并每个样本的 featureCounts 结果为矩阵（行=基因，列=样本，列顺序与样本表一致）。
mode:
  counts      -> counts.tsv + cpm.tsv
  length      -> tpm.tsv + fpkm.tsv（长度为 featureCounts 的 Length：基因非重叠外显子总长）
  counts_only -> counts.tsv
"""
import sys

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

FC_FIXED = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]

samples = list(snakemake.params.samples)
mode = snakemake.params.mode
files = list(snakemake.input)
assert len(files) == len(samples), "输入文件数与样本数不一致"

annot, counts = None, {}
for sample, fp in zip(samples, files):
    df = pd.read_csv(fp, sep="\t", comment="#", dtype={"Geneid": str})
    extra = [c for c in df.columns[len(FC_FIXED):-1]]  # --extraAttributes 列
    df = df.set_index("Geneid")
    if annot is None:
        annot = df[["Length", *extra]].rename(columns={"Length": "length"})
    elif not df.index.equals(annot.index):
        raise ValueError(f"{fp} 的基因顺序与其他样本不一致（是否使用了不同的 GTF？）")
    counts[sample] = df.iloc[:, -1]
    print(f"{sample}: {fp} assigned={int(df.iloc[:, -1].sum())}", file=sys.stderr)

counts = pd.DataFrame(counts)
counts.index.name = "gene_id"
annot.index.name = "gene_id"
# 注释列放前面: gene_id, gene_name, gene_biotype, ..., length
annot = annot[[c for c in annot.columns if c != "length"] + ["length"]]


def write(values, path, decimals=None):
    out = annot.join(values.round(decimals) if decimals is not None else values)
    out.to_csv(path, sep="\t")
    print(f"wrote {path} {out.shape}", file=sys.stderr)


if mode in ("counts", "counts_only"):
    write(counts, snakemake.output.counts)
if mode == "counts":
    cpm = counts / counts.sum().replace(0, float("nan")) * 1e6
    write(cpm.fillna(0), snakemake.output.cpm, 4)
if mode == "length":
    length = annot["length"]
    rpk = counts.div(length / 1e3, axis=0)
    tpm = rpk / rpk.sum().replace(0, float("nan")) * 1e6
    fpkm = counts.div(length, axis=0) * 1e9 / counts.sum().replace(0, float("nan"))
    write(tpm.fillna(0), snakemake.output.tpm, 4)
    write(fpkm.fillna(0), snakemake.output.fpkm, 4)
