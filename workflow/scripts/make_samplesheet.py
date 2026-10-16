#!/usr/bin/env python3
"""
扫描 fastq 目录生成 samples.tsv（sample / fq1 / fq2 + 可选的空白附加列）。

识别的文件名（.fastq.gz 或 .fq.gz）:
  X_1 / X_2,  X_R1 / X_R2,  X_R1_001 / X_R2_001,  X.R1 / X.R2   -> 双端 X
  其它                                                        -> 单端，样本名为去掉后缀的文件名
只有 _1 没有 _2（或反之）、或 _1/_2 的首条 read 名不一致（说明是两个单端重复）时，
按单端处理并保留完整文件名作为样本名。

用法: make_samplesheet.py -i <fq_dir> -o samples.tsv [--extra-cols strandedness,protocol,umi]
"""
import argparse
import gzip
import re
import sys
from pathlib import Path

SUFFIX = re.compile(r"\.(fastq|fq)\.gz$")
MATE = re.compile(r"^(?P<stem>.+?)[._](?P<r>R?[12])(?P<tail>_001)?$")


def first_read_id(path):
    with gzip.open(path, "rt") as f:
        head = f.readline().strip()
    rid = head[1:].split()[0] if head.startswith("@") else ""
    return re.sub(r"/[12]$", "", rid)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="fastq 所在目录（不递归）")
    ap.add_argument("-o", "--output", required=True, help="输出 samples.tsv")
    ap.add_argument("--extra-cols", default="", help="额外的空白列，逗号分隔")
    args = ap.parse_args()

    fq_dir = Path(args.input).resolve()
    files = sorted(p for p in fq_dir.iterdir() if p.is_file() and SUFFIX.search(p.name))
    if not files:
        sys.exit(f"错误: {fq_dir} 下没有 .fastq.gz / .fq.gz 文件")

    mates, singles = {}, []
    for p in files:
        stem = SUFFIX.sub("", p.name)
        m = MATE.match(stem)
        if m:
            mates.setdefault((m["stem"], m["tail"] or ""), {})[m["r"][-1]] = p
        else:
            singles.append((stem, p))

    rows = []
    for (stem, _), d in mates.items():
        if "1" in d and "2" in d and first_read_id(d["1"]) == first_read_id(d["2"]):
            rows.append((stem, d["1"], d["2"]))
        else:
            for p in d.values():
                print(f"提示: {p.name} 没有匹配的 mate，按单端处理", file=sys.stderr)
                singles.append((SUFFIX.sub("", p.name), p))
    rows += [(stem, p, "") for stem, p in singles]
    rows.sort(key=lambda r: r[0])

    names = [r[0] for r in rows]
    dup = sorted({n for n in names if names.count(n) > 1})
    if dup:
        sys.exit(f"错误: 样本名重复 {dup}（同一目录中可能同时存在 X.fq.gz 与 X.fastq.gz）")

    extra = [c for c in args.extra_cols.split(",") if c]
    with open(args.output, "w") as f:
        f.write("\t".join(["sample", "fq1", "fq2", *extra]) + "\n")
        for s, fq1, fq2 in rows:
            f.write("\t".join([s, str(fq1), str(fq2), *[""] * len(extra)]) + "\n")
    n_pe = sum(1 for r in rows if r[2])
    print(f"已生成 {args.output}: {len(rows)} 个样本（PE {n_pe}，SE {len(rows) - n_pe}）", file=sys.stderr)


if __name__ == "__main__":
    main()
