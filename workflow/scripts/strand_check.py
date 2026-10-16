"""
推断文库链方向：分别用 featureCounts -s 1 / -s 2 计数，比较 Assigned 数。
  frac_forward = assigned(-s 1) / (assigned(-s 1) + assigned(-s 2))
  >= threshold        -> forward
  <= 1 - threshold    -> reverse
  在 0.5 ± 0.1 之间    -> unstranded
  其他                -> undetermined
输出:
  {sample}.tsv : 推断明细
  {sample}.s   : 后续 featureCounts/StringTie 使用的 -s 值（auto 时用推断值，否则用填写值）
"""
import subprocess
import sys
import tempfile
from pathlib import Path

STRAND_TO_S = {"unstranded": 0, "forward": 1, "reverse": 2}

sys.stderr = open(snakemake.log[0], "w")

sample = snakemake.wildcards.sample
bam = snakemake.input.bam
gtf = snakemake.input.gtf
provided = snakemake.params.provided
threshold = float(snakemake.params.threshold)
pe_args = ["-p", "--countReadPairs"] if snakemake.params.pe else []


def assigned(s, tmpdir):
    out = Path(tmpdir) / f"s{s}.txt"
    cmd = [snakemake.params.featurecounts, "-T", str(snakemake.threads), "-s", str(s), *pe_args,
           "-t", "exon", "-g", "gene_id", "-a", gtf, "-o", str(out), bam]
    subprocess.run(cmd, check=True, stdout=sys.stderr, stderr=sys.stderr)
    for line in open(f"{out}.summary"):
        if line.startswith("Assigned"):
            return int(line.split("\t")[1])
    raise RuntimeError(f"{out}.summary 中没有 Assigned 行")


with tempfile.TemporaryDirectory(dir=Path(snakemake.output.tsv).parent) as tmp:
    a1, a2 = assigned(1, tmp), assigned(2, tmp)

total = a1 + a2
frac = a1 / total if total else float("nan")
if not total:
    inferred = "undetermined"
elif frac >= threshold:
    inferred = "forward"
elif frac <= 1 - threshold:
    inferred = "reverse"
elif abs(frac - 0.5) <= 0.1:
    inferred = "unstranded"
else:
    inferred = "undetermined"

if provided == "auto":
    used = inferred if inferred != "undetermined" else "unstranded"
    status = "auto" if inferred != "undetermined" else "auto_undetermined"
else:
    used = provided
    if inferred == "undetermined":
        status = "undetermined"
    else:
        status = "ok" if inferred == provided else "MISMATCH"

if status in ("MISMATCH", "auto_undetermined", "undetermined"):
    print(f"WARNING: {sample} provided={provided} inferred={inferred} "
          f"frac_forward={frac:.3f} -> 使用 {used}", file=sys.stderr)

with open(snakemake.output.tsv, "w") as f:
    f.write("sample\tprovided\tassigned_s1\tassigned_s2\tfrac_forward\tinferred\tused\tstatus\n")
    f.write(f"{sample}\t{provided}\t{a1}\t{a2}\t{frac:.4f}\t{inferred}\t{used}\t{status}\n")
with open(snakemake.output.s, "w") as f:
    f.write(f"{STRAND_TO_S[used]}\n")
