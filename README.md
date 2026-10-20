# snakemake_bulk_rna

bulk RNA-seq 本地分析流程（Snakemake ≥ 9.12）：

![workflow](docs/images/bulk_rna_workflow.svg)

<sub>实线框为主流程，虚线框为可选或按样本条件执行的步骤，绿色为 QC，★ 为交付结果（路径相对 `results/`）。PNG 版本：`docs/images/bulk_rna_workflow.png`。</sub>

<details>
<summary>文字版流程</summary>

```
fastq ─► (UMI extract) ─► trim_galore ─► hisat2 ─► (UMI dedup) ─► 链方向检查 ─► featureCounts ─► counts / CPM / TPM / FPKM
                 └─ FastQC (raw / trimmed)                     ├─► StringTie（可选）
                                                               ├─► repeat 计数（可选）
                                                               └─► bigwig
                                          汇总：summary.tsv + MultiQC 报告
```

</details>

支持按样本设置 **SE/PE**（自动识别）、**链方向**（可自动推断）、**全长 / 3' 端**、**有无 UMI**。

---

## 快速开始

```bash
P=/home/LAB/lab_ab_1/pipeline/snakemake_bulk_rna
cd /path/to/project

# 1. 生成项目配置模板（-c 指定的文件不存在时自动生成），修改其中的物种、文库类型等
bash $P/run.sh -c config.yaml

# 2. 扫描 fastq 目录，生成 <outdir>/samples.tsv，并预览任务
bash $P/run.sh -c config.yaml -i /path/to/rawdata -o rnaseq_out -n

# 3. 检查 rnaseq_out/samples.tsv 后，后台运行
nohup bash $P/run.sh -c config.yaml -o rnaseq_out -j 48 -y > run.log 2>&1 &
```

- 不需要激活任何环境，`run.sh` 直接调用 `~/.conda/envs/omics` 中的 snakemake（可用 `SMK_ENV=<环境路径>` 指定其它环境）。
- 流程代码只有一份，所有项目共用；**每个项目用 `-c` 指定自己的配置文件（必填）**。
  流程先读默认配置 `config/config.yaml`，再用 `-c` 的文件覆盖（嵌套的键逐层合并），所以项目配置可以只写与默认不同的参数。
  每次运行实际使用的（合并后的）配置会存档到 `<outdir>/logs/run_<时间>/config.yaml`。
  默认配置是所有项目的基础，一般不要改它；改了会影响所有项目的续跑。
- 结果、日志都写在 `-o` 指定的目录，多个项目可以同时跑。
- 中断或出错后，**重新运行同一条命令即可断点续跑**。

### run.sh 参数

| 参数 | 说明 |
|---|---|
| `-c FILE` | 项目配置文件（必需）。文件不存在时从默认配置生成模板后退出 |
| `-i DIR` | fastq 目录（仅首次运行、且 outdir 中还没有 samples.tsv 时使用） |
| `-o DIR` | 输出目录，默认 `./rnaseq_out` |
| `-s FILE` | 使用已写好的样本表 |
| `-j INT` | 最多同时使用的 CPU 核数，默认 32 |
| `-n` | 只预览（dry-run） |
| `-y` | 不询问，直接运行 |
| `-- ...` | 之后的参数原样传给 snakemake |

---

## 输入：samples.tsv

`run.sh -i` 自动生成，也可以手写（制表符分隔）：

```
sample   fq1                        fq2                        strandedness  protocol  umi
ctrl_1   /data/ctrl_1_R1.fq.gz      /data/ctrl_1_R2.fq.gz
ctrl_2   /data/ctrl_2.fastq.gz                                 reverse
q3_1     /data/q3_1.fastq.gz                                                 3prime    read
```

- `fq2` 为空即单端。识别的双端命名：`_1/_2`、`_R1/_R2`、`_R1_001/_R2_001`。
  若 `X_1`、`X_2` 的首条 read 名不同（说明是两个单端重复），会按两个单端样本处理。
- `strandedness` / `protocol` / `umi` **留空则使用配置中 `library:` 的默认值**，填写则覆盖该样本。
- 样本名只能包含字母、数字、`.` `_` `-`。

---

## 文库类型设置

### 1. 链方向 `strandedness`

| 值 | 含义 | 常见建库 |
|---|---|---|
| `reverse` | read1 与转录本反向（featureCounts `-s 2`） | Illumina TruSeq Stranded、NEBNext Ultra II Directional 等 dUTP 法；QuantSeq REV |
| `forward` | read1 与转录本同向（`-s 1`） | Lexogen QuantSeq FWD |
| `unstranded` | 无链信息（`-s 0`） | Smart-seq2、非链特异性 TruSeq / NEBNext Ultra |
| `auto` | 比对后自动推断（默认） | 不确定时使用 |

每个样本都会做一次链方向检查（分别以 `-s 1` 和 `-s 2` 计数）：

- 结果写在 `results/tables/strandedness.tsv`。
- `frac_forward ≥ 0.8` 判为 forward，`≤ 0.2` 判为 reverse，`0.4–0.6` 判为 unstranded。
- 填写值与推断值不一致时，`status` 列会标为 **MISMATCH**，流程仍然使用你填写的值。
- `auto` 时使用推断值；无法判断时按 unstranded 处理，`status` 为 `auto_undetermined`。

链方向会影响以下步骤：hisat2 `--rna-strandness`、featureCounts `-s`、StringTie `--fr/--rf`、分链 bigwig。

### 2. 全长 / 3' 端 `protocol`

| 值 | trim | 输出的表达矩阵 |
|---|---|---|
| `full_length` | trim_galore | counts、CPM、**TPM、FPKM**，可选 StringTie |
| `3prime` | trim_galore + 切除 read1 5' 端（`trim.3prime.clip_r1`）+ cutadapt 去 polyA | counts、CPM（3' 端测序按基因长度归一化没有意义，因此不输出 TPM/FPKM） |

### 3. UMI `umi`

| 值 | 处理方式 |
|---|---|
| `none` | 不去重（RNA-seq 常规做法） |
| `read` | UMI 在 read 序列中：比对前 `umi_tools extract`（按 `umi.pattern`），比对后 `umi_tools dedup` |
| `header` | UMI 已在 fastq 的 header 中：比对后 `umi_tools dedup`。UMI 取 header 中按 `umi.separator` 分隔的最后一个字段，位置由 `umi.header_field` 指定（见下表） |

| `umi.header_field` | UMI 位置 | 例子 |
|---|---|---|
| `name`（默认） | read 名末尾、空格之前 | bcl-convert：`@A00:...:2000:ACGTACGT 1:N:0:ATCACG` |
| `comment` | 空格后的注释末尾。比对软件会丢掉空格后的内容，所以比对前先把 read 名改写为 `@SRR123.1_ACGCAGGCCC` | SRA 下载的 DNBSEQ 数据：`@SRR123.1 V300...:2:CTCTCTAC:ACGCAGGCCC/1` |

流程启动时会抽查每个 `umi: header` 样本的前 1000 条 reads，取出的 UMI 不是等长的 ACGTN 序列时直接报错，避免因位置设错而静默地不去重。

有 UMI 时，后续定量使用去重后的 `results/03_align/{sample}.dedup.bam`。

### 常见建库的配置示例

以试剂盒说明书为准，以下仅作参考：

```yaml
# QuantSeq FWD（无 UMI）
library: {strandedness: forward, protocol: 3prime, umi: none}
trim: {3prime: {clip_r1: 12, polyA: true}}

# QuantSeq FWD + UMI Second Strand Synthesis Module（read1 开头 6nt UMI + 4nt TATA spacer）
library: {strandedness: forward, protocol: 3prime, umi: read}
umi: {pattern: NNNNNN}
trim: {3prime: {clip_r1: 4, polyA: true}}

# TruSeq Stranded mRNA，UMI 已由 bcl-convert 写入 read 名
library: {strandedness: reverse, protocol: full_length, umi: header}
umi: {header_field: name, separator: ":"}

# SRA 下载的 DNBSEQ 数据，UMI 在 header 空格后的注释末尾
library: {umi: header}
umi: {header_field: comment, separator: ":"}
```

---

## 配置 config.yaml 要点

| 配置项 | 说明 |
|---|---|
| `genome` | `mm10` / `GRCm39` / `GRCh38`，对应 `references` 中的 hisat2 index 与 Ensembl GTF |
| `gtf_filter.exclude_biotypes` | 定量前从 GTF 中去掉的基因类型，默认 `miRNA, snoRNA, snRNA`（与旧版 vM25 处理一致） |
| `threads` | 每类任务的线程数；总核数由 `run.sh -j` 控制，Snakemake 会自动安排并行 |
| `stringtie.enabled` | 额外输出 StringTie TPM/FPKM |
| `repeats.enabled` / `repeats.gtf` | repeat 计数（需自备 RepeatMasker GTF，染色体命名需与 hisat2 index 一致，mm10 index 为 `1, 2, …, MT`）。打开但 gtf 为空或不存在时，流程启动即报错 |
| `bigwig.split_strand` | 额外输出转录本正、负链 bigwig |
| `keep.trimmed_fastq` | 默认 false：trim 后的 fastq 用完即删 |
| `envs` | 软件所在的 conda 环境（`in202506`；deeptools、featureCounts 用 `omics`）。Python 脚本（汇总、合并矩阵等）直接用运行 snakemake 的环境，保证与 Snakemake 版本一致 |

流程启动时会检查 hisat2 index 与 GTF 的染色体命名是否一致，结果写在 `results/qc/reference_check.txt`。

---

## 输出

```
project/
├── samples.tsv                       # 本项目的样本表
├── results/
│   ├── tables/                        # ★ 最终交付
│   │   ├── counts_exon.tsv            # gene_id, gene_name, gene_biotype, length, 各样本 counts（只计 exon）
│   │   ├── cpm_exon.tsv
│   │   ├── counts_exon_intron.tsv     # exon + intron 计数（核 RNA 用）；length 为基因全长
│   │   ├── cpm_exon_intron.tsv
│   │   ├── tpm_exon.tsv / fpkm_exon.tsv   # 仅 full_length 样本，按 exon 长度
│   │   ├── stringtie_tpm.tsv / stringtie_fpkm.tsv   # 可选
│   │   ├── repeats_counts.tsv         # 可选
│   │   ├── strandedness.tsv           # 链方向推断
│   │   └── summary.tsv                # 每个样本的 reads 数、比对率、assigned 率等
│   ├── 00_raw/                        # 原始 fastq 的软链接（统一命名为 {sample}_R1/R2）
│   ├── 02_trim/reports/               # trim_galore 报告
│   ├── 03_align/                      # {sample}.bam（+ .dedup.bam）
│   ├── 04_count/                      # 每个样本的 featureCounts 结果、strand/
│   ├── 05_stringtie/  06_repeats/  07_bigwig/
│   └── qc/
│       ├── multiqc/multiqc_report.html  # ★ 所有 QC 汇总
│       ├── fastqc_raw/  fastqc_trimmed/
│       ├── hisat2/{sample}.summary.txt  # 比对统计
│       ├── umi/                         # umi_tools extract / dedup 统计（仅 UMI 样本）
│       └── reference_check.txt          # index 与 GTF 染色体命名检查
└── logs/
    └── run_<时间>/                    # 每执行一次 snakemake 新建一个；续跑时只含这次实际跑了的任务
        ├── <步骤>/<样本>.log          # 出错时先看这里
        ├── benchmarks.tsv             # 本次所有任务的耗时与内存（一张表）
        ├── snakemake.log              # snakemake 完整输出
        ├── config.yaml、samples.tsv   # 本次使用的配置与样本表
        ├── pipeline.txt               # 流程路径、commit 版本、运行命令
        └── software_versions.tsv      # 各软件的版本与实际路径
```

**核 RNA（nuclear RNA-seq）**：比对不变，所有样本同时输出两套计数，按数据类型选用：

- `counts_exon.tsv` / `cpm_exon.tsv`：featureCounts `-t exon -g gene_id`，全细胞 RNA 使用。
- `counts_exon_intron.tsv` / `cpm_exon_intron.tsv`：参考 zUMIs，exon 优先。先按 exon 计数，未分配给 exon 的 reads
  再按基因全长（`-t gene`）计数并加到 exon 计数上，一条 read 只算一次。核 RNA 含大量未剪接的 pre-mRNA，使用这一套。
- 两个基因重叠区域（嵌套基因、反义基因）的 reads 会被判为 ambiguous，不计入任何基因，与 exon 计数的规则一致。
- 内含子计数依赖链方向；非链特异的数据中，重叠基因的内含子 reads 损失会更多。
- TPM / FPKM 只按 exon 计算。

**summary.tsv 的各列**（reads 按单条 read 计，PE 为 R1 + R2；所有 rate 都以 `raw_reads` 为分母，单位 %；reads 数带千位分隔符）：

| 列 | 含义 |
|---|---|
| `raw_reads` | trim 前的 reads |
| `clean_reads` | 进入 hisat2 的 reads |
| `mapped_reads` | 比对上的 reads（clean_reads × hisat2 overall alignment rate） |
| `umi_dedup_reads` | UMI 去重后的 reads（仅 UMI 样本）。按 read 名计数，多重比对的 read 只算一次。MultiQC 的 UMI-tools 部分是 umi_tools 自己按比对记录计的，数值会偏大，以本表为准 |
| `clean_rate (%)` | clean_reads / raw_reads |
| `mapping_rate (%)` | mapped_reads / raw_reads（与 hisat2 日志中以 clean 为分母的比对率不同） |
| `umi_dedup_rate (%)` / `umi_duplicates_rate (%)` | 去重后 reads / raw_reads；(去重前比对上的 reads − 去重后) / raw_reads（仅 UMI 样本） |
| `assigned_rate (%)` | featureCounts 计入基因的 reads / raw_reads（PE 的片段数 × 2） |
| `intron_rate (%)` | 内含子 reads / raw_reads（见下方"核 RNA"）。全细胞 RNA 通常较低，核 RNA 明显更高 |
| `strand_*` | 链方向：填写值、推断值、实际使用值、状态 |

---

## 常用操作

```bash
# 断点续跑：重新运行同一条命令
bash $P/run.sh -c config.yaml -o rnaseq_out -y

# 只跑到某一步（所有样本）
bash $P/run.sh -c config.yaml -o rnaseq_out -y -- --until featurecounts

# 强制重跑某一步（及其下游），例如改了 featureCounts 参数
bash $P/run.sh -c config.yaml -o rnaseq_out -y -- --forcerun featurecounts

# 追加样本：在 samples.tsv 中加行后直接运行，已完成的样本不会重跑

# 上次运行被强行中断，提示目录被锁
bash $P/run.sh -c config.yaml -o rnaseq_out -y -- --unlock

# 查看每一步为什么要（重）跑
bash $P/run.sh -c config.yaml -o rnaseq_out -n -- --summary
```

修改 config 中的参数或更新流程代码后，Snakemake 会自动重跑受影响的步骤。

---

## 与旧版的主要区别

- 不再提交 PBS，改为本地运行：每条规则声明线程数，由 Snakemake 统一调度。
- 不再通过 `--config reads=<JSON>` 传参，改用可编辑的 samples.tsv；支持 `.fq.gz`。
- SE/PE 规则合并，文件名不再带 `_pe` / `_se`。
- 定量改为按 `gene_id` 计数，矩阵附带 `gene_name` 等注释列。cufflinks 被 featureCounts 计算的 TPM/FPKM 取代，另有可选的 StringTie。
- 修复的问题：
  - repeat 合并结果实际上是基因 count；
  - summary 读不到 snakemake 参数；
  - human/mm9 的 GTF 配置结构错误；
  - 线程数写死。
- 新增：链方向推断、3' 端测序与 UMI 支持、统一的 MultiQC 报告、每次运行的记录 `logs/run_<时间>/`。
