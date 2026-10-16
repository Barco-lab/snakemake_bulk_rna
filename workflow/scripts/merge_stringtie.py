"""合并 StringTie gene_abund.tab 为 TPM / FPKM 矩阵。同一 Gene ID 出现多行时求和。"""
import sys

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

names, tpm, fpkm = {}, {}, {}
for sample, fp in zip(snakemake.params.samples, snakemake.input):
    df = pd.read_csv(fp, sep="\t", dtype={"Gene ID": str})
    names.update(df.drop_duplicates("Gene ID").set_index("Gene ID")["Gene Name"].to_dict())
    g = df.groupby("Gene ID")[["TPM", "FPKM"]].sum()
    tpm[sample], fpkm[sample] = g["TPM"], g["FPKM"]
    print(f"{sample}: {fp} genes={len(g)}", file=sys.stderr)

for values, path in ((tpm, snakemake.output.tpm), (fpkm, snakemake.output.fpkm)):
    out = pd.DataFrame(values).fillna(0).round(4)
    out.insert(0, "gene_name", pd.Series(names).reindex(out.index))
    out.index.name = "gene_id"
    out.sort_index().to_csv(path, sep="\t")
