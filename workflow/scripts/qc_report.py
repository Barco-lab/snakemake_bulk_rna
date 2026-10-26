"""
单文件 HTML QC 报告：概览、流程与参数、QC 汇总表（中英文可切换）。
所有样式、脚本和流程图都内嵌，不依赖网络。各列说明、参考值和检查规则来自 workflow/report/qc_metrics.yaml：
  checks  有公开阈值的检查（warn 橙 / fail 红）
  outlier 与本批样本比较，|值 - 中位数| > 3 x MAD（x 1.4826）标为橙色
"""
import html
import re
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml

sys.stderr = open(snakemake.log[0], "w")

p = snakemake.params
df = pd.read_csv(snakemake.input.summary, sep="\t")
strand = pd.read_csv(snakemake.input.strand, sep="\t").set_index("sample")
meta = yaml.safe_load(Path(snakemake.input.metrics).read_text())
METRICS, SOURCES = meta["metrics"], meta.get("sources", {})


def base(col):
    """summary.tsv 列名去掉括号部分，如 'mapping_rate (% of clean)' -> 'mapping_rate'。"""
    return col.split(" (")[0]


def esc(x):
    return html.escape(str(x))


def rich(text):
    """转义后把 `...` 显示为代码。"""
    return re.sub(r"`([^`]+)`", r"<code>\1</code>", esc(text))


def L(zh, en=None):
    """中英文各一份，按当前语言显示其一。"""
    en = zh if en is None else en
    if zh == en:
        return rich(zh)
    return f'<span data-l="zh">{rich(zh)}</span><span data-l="en">{rich(en)}</span>'


def T(zh, en):
    """悬停说明（title），切换语言时由脚本替换。"""
    return f' title="{esc(zh)}" data-tz="{esc(zh)}" data-te="{esc(en)}"'


def zh_en(x):
    return (x["zh"], x["en"]) if isinstance(x, dict) else (x, x)


def is_num(col):
    return pd.api.types.is_numeric_dtype(df[col])


def fmt(col, v):
    """reads 等计数加千位分隔符；其余数值按 qc_metrics.yaml 的 digits（默认 2）位小数。"""
    if pd.isna(v):
        return ""
    if base(col).endswith("_reads") or pd.api.types.is_integer_dtype(df[col]):
        return f"{int(v):,}"
    if is_num(col):
        return f"{v:.{METRICS.get(base(col), {}).get('digits', 2)}f}"
    return str(v)


cols = {base(c): c for c in df.columns}
col = lambda key: df[cols[key]] if key in cols else None
n = len(df)

# ------------------------------------------------------------
# 检查：flags[(行号, 列名)] = [(level, 中文, 英文)]
# ------------------------------------------------------------
flags = {}


def flag(i, c, level, zh, en):
    flags.setdefault((i, c), []).append((level, zh, en))


for i, row in df.iterrows():
    env = {base(c): (None if pd.isna(v) else v) for c, v in row.items()}
    for key, m in METRICS.items():
        if key not in cols:
            continue
        for chk in m.get("checks", []):
            if any(env.get(k) != v for k, v in chk.get("when", {}).items()):
                continue
            try:
                hit = eval(chk["expr"], {}, env)
            except TypeError:  # 缺失值
                continue
            if hit:
                flag(i, cols[key], chk["level"], *zh_en(chk["message"]))

outliers = []  # (样本, 列, 高/低, 中位数)
for key, m in METRICS.items():
    if not m.get("outlier") or key not in cols:
        continue
    sel = df
    for k, v in m.get("outlier_when", {}).items():  # 只在部分样本之间比较
        sel = sel[sel[cols[k]] == v]
    s = sel[cols[key]].dropna().astype(float)
    if len(s) < 4:
        continue
    med = s.median()
    mad = 1.4826 * (s - med).abs().median()
    if mad == 0:
        continue
    for i, v in s.items():
        if abs(v - med) > 3 * mad:
            hi = v > med
            outliers.append((df.loc[i, "sample"], key, hi, fmt(cols[key], med)))
            mt = fmt(cols[key], med)
            flag(i, cols[key], "warn",
                 f"明显{'高' if hi else '低'}于本批中位数 {mt}（|差| > 3 × MAD）",
                 f"{'above' if hi else 'below'} the batch median {mt} (|diff| > 3 × MAD)")

flagged = sorted({i for i, _ in flags})

# 指定了链方向时，推断值与实际使用值不显示（推断结果体现在 strand_status 中）
hidden = set()
if "strand_provided" in cols and (col("strand_provided") != "auto").all():
    hidden = {cols[k] for k in ("strand_inferred", "strand_used") if k in cols}
shown = [c for c in df.columns if c not in hidden]


# ------------------------------------------------------------
# 概览
# ------------------------------------------------------------
def med(key, scale=1, digits=1, unit=""):
    s = col(key)
    if s is None or s.dropna().empty:
        return None
    return f"{s.median() / scale:.{digits}f}{unit}"


ok_status = ("ok", "auto")
strand_ok = int(col("strand_status").isin(ok_status).sum()) if "strand_status" in cols else None
tiles = [
    (("样本", "Samples"), str(n), ("", "")),
    (("raw reads 中位数", "Median raw reads"), med("raw_reads", 1e6, 1, "M"), ("", "")),
    (("mapping rate 中位数", "Median mapping rate"), med("mapping_rate", 1, 1, "%"), ("以 clean reads 为分母", "of clean reads")),
    (("assigned rate 中位数", "Median assigned rate"), med("assigned_rate", 1, 1, "%"), ("以进入 featureCounts 的 reads 为分母", "of reads counted")),
    (("UMI 去重后 reads 中位数", "Median UMI-dedup reads"), med("umi_dedup_reads", 1e6, 1, "M"), ("", "")),
    (("链方向一致", "Strandedness consistent"), f"{strand_ok} / {n}" if strand_ok is not None else None, ("", "")),
    (("需关注的样本", "Samples to review"), str(len(flagged)), ("见下方列表", "listed below")),
]
tiles_html = "".join(
    f'<div class="tile{" tile-alert" if label[0] == "需关注的样本" and flagged else ""}">'
    f'<div class="tile-label">{L(*label)}</div><div class="tile-value">{esc(value)}</div>'
    f'<div class="tile-note">{L(*note)}</div></div>'
    for label, value, note in tiles
    if value is not None
)

if flagged:
    items = []
    for i in flagged:
        chips = "".join(
            f'<span class="chip chip-{lvl}">{esc(base(c))}：{L(zh, en)}</span>'
            for (j, c), lst in flags.items()
            if j == i
            for lvl, zh, en in lst
        )
        items.append(f'<tr><td class="sample">{esc(df.loc[i, "sample"])}</td><td>{chips}</td></tr>')
    attention = f'<table class="plain"><tbody>{"".join(items)}</tbody></table>'
else:
    attention = f'<p class="ok-note">{L("所有样本均未触发检查，也没有离群值。", "No sample triggered a check or was an outlier.")}</p>'


# ------------------------------------------------------------
# 要点：根据数据自动生成，放在 QC 汇总表下方
# ------------------------------------------------------------
def rng(s, scale=1, digits=1, unit=""):
    s = s.dropna() / scale
    return f"{s.median():.{digits}f}{unit}", f"{s.min():.{digits}f}–{s.max():.{digits}f}{unit}"


def names(lst):
    return ", ".join(lst)


points = []  # (中文, 英文)
prot = col("protocol")
raw_m, raw_r = rng(col("raw_reads"), 1e6, 1, "M")
zh, en = f"测序深度：raw reads 中位数 {raw_m}（{raw_r}）。", f"Depth: median raw reads {raw_m} ({raw_r})."
three = df[prot == "3prime"]
if len(three):
    low = three.loc[three[cols["raw_reads"]] < 3e6, "sample"].tolist()
    zh += f" 3' 端测序样本{'均达到' if not low else '中有 ' + str(len(low)) + ' 个低于'} Lexogen 建议的 3M reads" + (f"（{names(low)}）" if low else "") + "。"
    en += f" {'All' if not low else str(len(low))} 3'-end samples {'meet' if not low else 'fall below'} Lexogen's 3M-read recommendation" + (f" ({names(low)})" if low else "") + "."
full = df[prot == "full_length"]
if len(full):
    low = full.loc[full[cols["mapped_reads"]] < 30e6, "sample"].tolist()
    zh += f" full_length 样本{'均达到' if not low else '中有 ' + str(len(low)) + ' 个低于'} ENCODE 30M aligned reads" + (f"（{names(low)}）" if low else "") + "。"
    en += f" {'All' if not low else str(len(low))} full-length samples {'meet' if not low else 'fall below'} the ENCODE 30M aligned-read standard" + (f" ({names(low)})" if low else "") + "."
if "umi_dedup_reads" in cols and col("umi_dedup_reads").notna().any():
    d_m, d_r = rng(col("umi_dedup_reads"), 1e6, 1, "M")
    u_m, u_r = rng(col("umi_duplicates_rate"))
    zh += f" UMI 去重后中位数 {d_m}（{d_r}），UMI 重复率中位数 {u_m}%（{u_r}%），独立分子数应按去重后的 reads 评估。"
    en += f" After UMI deduplication the median is {d_m} ({d_r}) with a median UMI duplicate rate of {u_m}% ({u_r}%); judge library depth by deduplicated reads."
points.append((zh, en))

h_m, h_r = rng(col("mapping_rate"))
low = df.loc[col("mapping_rate") < 70, "sample"].tolist()
points.append((
    f"比对：mapping rate（以 clean 为分母）中位数 {h_m}%（{h_r}%），"
    + ("均在 Conesa 2016 给出的常规范围（≥ 70%）内。" if not low else f"{len(low)} 个样本低于 70%（{names(low)}），需检查污染或 rRNA。"),
    f"Alignment: median mapping rate (of clean reads) {h_m}% ({h_r}%), "
    + ("all within the typical range (≥ 70%, Conesa 2016)." if not low else f"{len(low)} below 70% ({names(low)}); check contamination or rRNA."),
))

a_m, a_r = rng(col("assigned_rate"))
i_m, i_r = rng(col("intron_rate"))
points.append((
    f"定量：进入 featureCounts 的 reads 中，assigned rate 中位数 {a_m}%（{a_r}%），intron rate 中位数 {i_m}%（{i_r}%）。"
    "polyA / 3' 端文库的内含子比例通常较低，核 RNA / 总 RNA 更高；本应较低时偏高，提示基因组 DNA 污染。",
    f"Quantification: of the reads entering featureCounts, median assigned rate {a_m}% ({a_r}%) and intron rate {i_m}% ({i_r}%). "
    "Intronic fractions are usually low for polyA / 3'-end libraries and higher for nuclear / total RNA; unexpectedly high values suggest gDNA contamination.",
))

if "strand_status" in cols:
    st = col("strand_status")
    und = df.loc[st.isin(["undetermined", "auto_undetermined"]), "sample"].tolist()
    mis = df.loc[st == "MISMATCH", "sample"].tolist()
    zh = f"链方向：{strand_ok}/{n} 个样本推断结果与填写值一致。"
    en = f"Strandedness: {strand_ok}/{n} samples agree with the given value."
    if und:
        fr = strand.loc[und, "frac_forward"]
        zh += f" {len(und)} 个样本无法判断（forward 比例 {fr.min():.2f}–{fr.max():.2f}，未达阈值）：{names(und)}。"
        en += f" {len(und)} undetermined (forward fraction {fr.min():.2f}–{fr.max():.2f}, below threshold): {names(und)}."
        if hidden:  # 链方向已指定
            zh += " 由于已指定链方向，定量仍按填写值进行；如比例接近 0.5，需确认文库是否确为链特异。"
            en += " Quantification still uses the given strandedness; if the fraction is close to 0.5, confirm the library is stranded."
    if mis:
        zh += f" {len(mis)} 个样本与填写值相反（{names(mis)}），请检查样本表。"
        en += f" {len(mis)} contradict the given value ({names(mis)}); check the sample sheet."
    points.append((zh, en))

if outliers:
    zh = "离群值（与本批中位数相差超过 3 × MAD）：" + "；".join(
        f"{s} 的 {k}{'偏高' if hi else '偏低'}（中位数 {m}）" for s, k, hi, m in outliers) + "。"
    en = "Outliers (> 3 × MAD from the batch median): " + "; ".join(
        f"{s} {k} {'high' if hi else 'low'} (median {m})" for s, k, hi, m in outliers) + "."
    points.append((zh, en))

FAIL_NOTE = ("建议重点检查或剔除。", "Review or exclude them.")
fails = {}  # 列名 -> 低于公开标准的样本
for (i, c), lst in flags.items():
    if any(l == "fail" for l, *_ in lst):
        fails.setdefault(base(c), []).append(df.loc[i, "sample"])
if fails:
    zh = "总体：低于公开标准的样本——" + "；".join(f"{k}：{names(v)}" for k, v in fails.items()) + "。" + FAIL_NOTE[0]
    en = "Overall: samples below published standards — " + "; ".join(f"{k}: {names(v)}" for k, v in fails.items()) + ". " + FAIL_NOTE[1]
else:
    zh = "总体：没有样本低于公开标准；" + ("标为橙色的样本建议结合实验信息复核。" if flagged else "整体质量良好。")
    en = "Overall: no sample falls below published standards; " + ("review the samples marked in orange with the experimental context." if flagged else "overall quality is good.")
points.append((zh, en))
conclusion = "".join(f"<li>{L(zh, en)}</li>" for zh, en in points)
sources = " · ".join(
    f'<a href="{esc(s["url"])}" target="_blank" rel="noopener" title="{esc(s["name"])}">{esc(k)}</a>'
    for k, s in SOURCES.items()
)

# ------------------------------------------------------------
# 流程与参数
# ------------------------------------------------------------
svg = re.sub(r"^<\?xml[^>]*>\s*", "", Path(snakemake.input.flowchart).read_text())
settings_html = "".join(f"<tr><th>{L(kz, ke)}</th><td>{L(vz, ve)}</td></tr>" for kz, ke, vz, ve in p.settings)


# ------------------------------------------------------------
# QC 汇总表
# ------------------------------------------------------------
def head_tip(c):
    m = METRICS.get(base(c), {})
    dz, de = zh_en(m.get("desc", ""))
    if m.get("ref"):
        rz, re_ = zh_en(m["ref"])
        dz += f"\n参考值：{rz}"
        de += f"\nReference: {re_}"
    if m.get("sources"):
        dz += "\n出处：" + ", ".join(m["sources"])
        de += "\nSource: " + ", ".join(m["sources"])
    return T(dz, de) if dz else ""


head = "".join(
    f'<th{head_tip(c)}{" class=num" if is_num(c) else ""}>{esc(c)}<span class="sort"></span></th>' for c in shown
)
body = []
for i, row in df.iterrows():
    cells = []
    for c in shown:
        v = row[c]
        cls = (["num"] if is_num(c) else []) + (["sample"] if c == "sample" else [])
        tip = ""
        if (i, c) in flags:
            lst = flags[(i, c)]
            cls.append("fail" if any(l == "fail" for l, *_ in lst) else "warn")
            tip = T("；".join(z for _, z, _ in lst), "; ".join(e for *_, e in lst))
        data = "" if pd.isna(v) else str(v)
        cells.append(f'<td class="{" ".join(cls)}" data-v="{esc(data)}"{tip}>{esc(fmt(c, v))}</td>')
    body.append("<tr>" + "".join(cells) + "</tr>")
table_html = f'<table class="qc" id="qc" data-key="{esc(p.pipeline)}"><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table>'

# ------------------------------------------------------------
# 页面
# ------------------------------------------------------------
title_zh, title_en = p.title
# 项目：输出目录（snakemake --directory，即当前目录）的上两级，如 Patri_aging/out_rna
OUTDIR = Path.cwd()
PROJECT = "/".join(OUTDIR.parts[-2:])
CSS = """
:root {
  --bg: #f5f5f7; --card: #ffffff; --text: #1d1d1f; --muted: #6e6e73; --line: #d2d2d7; --soft: #f5f5f7;
  --accent: #0071e3; --green: #248a3d; --orange: #c93400; --red: #d70015;
  --orange-bg: #fff4e5; --red-bg: #ffebee; --shadow: 0 2px 12px rgba(0,0,0,.06); --nav: rgba(255,255,255,.72);
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #000000; --card: #1c1c1e; --text: #f5f5f7; --muted: #98989d; --line: #38383a; --soft: #2c2c2e;
    --accent: #2997ff; --green: #30d158; --orange: #ff9f0a; --red: #ff453a;
    --orange-bg: rgba(255,159,10,.16); --red-bg: rgba(255,69,58,.18); --shadow: none; --nav: rgba(28,28,30,.72);
  }
}
* { box-sizing: border-box; }
html { scroll-behavior: smooth; scroll-padding-top: 24px; }
body {
  margin: 0; background: var(--bg); color: var(--text);
  font-family: -apple-system, BlinkMacSystemFont, "SF Pro Text", "PingFang SC", "Helvetica Neue", "Noto Sans CJK SC", "Microsoft YaHei", sans-serif;
  font-size: 14px; line-height: 1.55; -webkit-font-smoothing: antialiased;
}
body.en [data-l="zh"], body:not(.en) [data-l="en"] { display: none; }
a { color: var(--accent); text-decoration: none; }
a:hover { text-decoration: underline; }
code { font-family: "SF Mono", ui-monospace, Menlo, Consolas, monospace; font-size: 12.5px; background: var(--soft); padding: 1px 6px; border-radius: 6px; word-break: break-all; }
nav {
  position: fixed; inset: 0 auto 0 0; width: 232px; padding: 28px 16px; background: var(--nav);
  backdrop-filter: saturate(180%) blur(20px); -webkit-backdrop-filter: saturate(180%) blur(20px); border-right: 1px solid var(--line);
}
.brand { font-size: 17px; font-weight: 600; letter-spacing: -.01em; padding: 0 12px; }
.brand-sub { color: var(--muted); font-size: 12px; padding: 2px 12px 16px; }
.seg { display: inline-flex; background: var(--soft); border-radius: 9px; padding: 2px; margin: 0 12px 16px; }
.seg button { border: 0; background: none; color: var(--text); font: inherit; font-size: 12.5px; padding: 4px 14px; border-radius: 7px; cursor: pointer; }
.seg button.on { background: var(--card); box-shadow: 0 1px 3px rgba(0,0,0,.12); font-weight: 600; }
nav a.item { display: block; padding: 7px 12px; margin: 2px 0; border-radius: 8px; color: var(--text); }
nav a.item:hover { background: var(--soft); text-decoration: none; }
nav a.item.active { background: var(--accent); color: #fff; }
nav .ext { margin-top: 18px; padding-top: 14px; border-top: 1px solid var(--line); }
nav .ext a { display: block; padding: 6px 12px; font-size: 13px; }
main { margin-left: 232px; padding: 40px 48px 80px; max-width: 1480px; }
header.top h1 { font-size: 32px; font-weight: 700; letter-spacing: -.02em; margin: 0; }
header.top .project { margin-top: 8px; font-size: 15px; font-weight: 600; }
header.top .project .path { color: var(--muted); font-weight: 400; font-size: 13px; margin-left: 8px; overflow-wrap: anywhere; }
header.top .meta { color: var(--muted); margin-top: 6px; font-size: 13px; }
section { margin-top: 44px; }
section > h2 { font-size: 24px; font-weight: 600; letter-spacing: -.015em; margin: 0 0 16px; }
h3 { font-size: 17px; font-weight: 600; margin: 0 0 12px; }
.card { background: var(--card); border-radius: 18px; padding: 24px; box-shadow: var(--shadow); margin-bottom: 20px; }
.tiles { display: grid; grid-template-columns: repeat(auto-fill, minmax(180px, 1fr)); gap: 14px; margin-bottom: 20px; }
.tile { background: var(--card); border-radius: 16px; padding: 16px 18px; box-shadow: var(--shadow); }
.tile-label { color: var(--muted); font-size: 12.5px; }
.tile-value { font-size: 28px; font-weight: 600; letter-spacing: -.02em; margin-top: 4px; font-variant-numeric: tabular-nums; }
.tile-note { color: var(--muted); font-size: 11.5px; min-height: 1.2em; }
.tile-alert .tile-value { color: var(--orange); }
.chip { display: inline-block; font-size: 12.5px; padding: 3px 10px; border-radius: 999px; margin: 2px 6px 2px 0; }
.chip-warn { background: var(--orange-bg); color: var(--orange); }
.chip-fail { background: var(--red-bg); color: var(--red); }
.ok-note { color: var(--green); margin: 0; }
table { border-collapse: collapse; width: 100%; }
table.plain td { padding: 8px 0; border-bottom: 1px solid var(--line); vertical-align: top; }
table.plain tr:last-child td { border-bottom: 0; }
table.plain td.sample { width: 220px; font-weight: 500; padding-right: 16px; }
table.kv th { text-align: left; font-weight: 500; color: var(--muted); width: 200px; padding: 10px 16px 10px 0; vertical-align: top; }
table.kv td { padding: 10px 0; }
table.kv tr + tr th, table.kv tr + tr td { border-top: 1px solid var(--line); }
.flow { background: #fff; border-radius: 12px; padding: 8px; max-width: 760px; margin: 0 auto; }
.flow svg { width: 100%; height: auto; display: block; }
.muted { color: var(--muted); font-size: 12.5px; }
.legend span.sw { display: inline-block; width: 10px; height: 10px; border-radius: 3px; margin: 0 4px 0 12px; vertical-align: -1px; }
.bar { display: flex; align-items: center; justify-content: space-between; gap: 12px; margin-bottom: 12px; }
.btn { border: 0; background: var(--accent); color: #fff; font: inherit; font-size: 13px; padding: 6px 14px; border-radius: 999px; cursor: pointer; white-space: nowrap; }
.btn:hover { filter: brightness(1.08); }
.tools { display: flex; gap: 8px; position: relative; }
.btn.ghost { background: var(--soft); color: var(--text); border: 1px solid var(--line); }
.colpick { display: none; position: absolute; right: 0; top: calc(100% + 6px); z-index: 10; background: var(--card);
  border: 1px solid var(--line); border-radius: 12px; box-shadow: var(--shadow); padding: 10px 14px;
  max-height: 60vh; overflow: auto; min-width: 220px; font-size: 13px; }
.colpick.open { display: block; }
.colpick label { display: block; white-space: nowrap; padding: 2px 0; cursor: pointer; }
.colpick .all { color: var(--accent); cursor: pointer; margin-bottom: 6px; display: inline-block; }
.scroll { overflow: auto; max-height: 78vh; border-radius: 12px; border: 1px solid var(--line); }
table.qc { font-size: 13px; font-variant-numeric: tabular-nums; }
table.qc th {
  position: sticky; top: 0; z-index: 2; background: var(--soft); text-align: left; font-weight: 600; white-space: nowrap;
  padding: 10px 12px; border-bottom: 1px solid var(--line); cursor: pointer; user-select: none;
}
table.qc th:hover { color: var(--accent); }
table.qc td { padding: 8px 12px; border-bottom: 1px solid var(--line); white-space: nowrap; background: var(--card); }
table.qc tr:last-child td { border-bottom: 0; }
table.qc tbody tr:hover td { background: var(--soft); }
table.qc .num { text-align: right; }
table.qc th:first-child, table.qc td:first-child { position: sticky; left: 0; z-index: 1; }
table.qc th:first-child { z-index: 3; }
table.qc td.sample { font-weight: 500; }
table.qc td.warn { background: var(--orange-bg); color: var(--orange); font-weight: 600; }
table.qc td.fail { background: var(--red-bg); color: var(--red); font-weight: 600; }
.sort::after { content: ""; margin-left: 4px; color: var(--accent); }
th[data-dir="asc"] .sort::after { content: "▲"; font-size: 9px; }
th[data-dir="desc"] .sort::after { content: "▼"; font-size: 9px; }
ul.points { margin: 0; padding-left: 20px; }
ul.points li { margin: 8px 0; }
ul.points li:last-child { font-weight: 600; }
.digest { border-left: 4px solid var(--accent); }
@media (max-width: 860px) {
  nav { position: static; width: auto; border-right: 0; border-bottom: 1px solid var(--line); padding: 16px; }
  nav a.item { display: inline-block; }
  nav .ext { display: none; }
  main { margin-left: 0; padding: 24px 16px 60px; }
  table.kv th { width: 120px; }
}
"""
JS = """
// 中英文切换（记住上次的选择）
function setLang(l) {
  document.body.classList.toggle("en", l === "en");
  document.documentElement.lang = l === "en" ? "en" : "zh-CN";
  document.querySelectorAll("[data-tz]").forEach(e => e.title = l === "en" ? e.dataset.te : e.dataset.tz);
  document.querySelectorAll(".seg button").forEach(b => b.classList.toggle("on", b.dataset.lang === l));
  try { localStorage.setItem("qc-lang", l); } catch (e) {}
}
document.querySelectorAll(".seg button").forEach(b => b.addEventListener("click", () => setLang(b.dataset.lang)));
let saved = null;
try { saved = localStorage.getItem("qc-lang"); } catch (e) {}
setLang(saved || "zh");

// 表头排序（按 data-v，空值始终在最后）
document.querySelectorAll("table.qc th").forEach((th, idx) => th.addEventListener("click", () => {
  const table = th.closest("table"), tbody = table.tBodies[0];
  const dir = th.dataset.dir === "asc" ? "desc" : "asc";
  table.querySelectorAll("th").forEach(h => delete h.dataset.dir);
  th.dataset.dir = dir;
  const val = r => { const v = r.cells[idx].dataset.v; const x = Number(v); return v === "" ? null : (isNaN(x) ? v : x); };
  [...tbody.rows].sort((a, b) => {
    const x = val(a), y = val(b);
    if (x === null) return 1; if (y === null) return -1;
    return (x > y ? 1 : x < y ? -1 : 0) * (dir === "asc" ? 1 : -1);
  }).forEach(r => tbody.appendChild(r));
}));

// 显示 / 隐藏列（sample 列始终显示；选择按流程记在浏览器中）
(() => {
  const t = document.getElementById("qc"), box = document.getElementById("colpick");
  const key = "qc-cols-" + t.dataset.key;
  const names = [...t.tHead.rows[0].cells].map(c => c.firstChild.textContent);
  let hidden = new Set();
  try { hidden = new Set(JSON.parse(localStorage.getItem(key) || "[]")); } catch (e) {}
  const apply = () => {
    [...t.rows].forEach(r => [...r.cells].forEach((c, i) => { c.hidden = i > 0 && hidden.has(names[i]); }));
    box.querySelectorAll("input").forEach(x => { x.checked = !hidden.has(x.value); });
    try { localStorage.setItem(key, JSON.stringify([...hidden])); } catch (e) {}
  };
  names.slice(1).forEach(n => {
    const lab = document.createElement("label"), x = document.createElement("input");
    x.type = "checkbox"; x.value = n;
    x.addEventListener("change", () => { x.checked ? hidden.delete(n) : hidden.add(n); apply(); });
    lab.append(x, " " + n); box.appendChild(lab);
  });
  box.querySelector(".all").addEventListener("click", () => { hidden.clear(); apply(); });
  document.getElementById("cols").addEventListener("click", e => { e.stopPropagation(); box.classList.toggle("open"); });
  box.addEventListener("click", e => e.stopPropagation());
  document.addEventListener("click", () => box.classList.remove("open"));
  apply();
})();

// 复制表格：制表符分隔、纯数字，可直接粘贴到 Excel（只复制显示的列）
document.getElementById("copy").addEventListener("click", async e => {
  const t = document.getElementById("qc");
  const rows = [[...t.tHead.rows[0].cells].filter(c => !c.hidden).map(c => c.textContent)];
  [...t.tBodies[0].rows].forEach(r => rows.push([...r.cells].filter(c => !c.hidden).map(c => c.dataset.v)));
  const text = rows.map(r => r.join("\\t")).join("\\n");
  try { await navigator.clipboard.writeText(text); }
  catch (err) {
    const ta = document.createElement("textarea"); ta.value = text; document.body.appendChild(ta);
    ta.select(); document.execCommand("copy"); ta.remove();
  }
  const b = e.currentTarget; b.classList.add("done");
  setTimeout(() => b.classList.remove("done"), 1500);
});

// 左侧导航随滚动高亮
const links = [...document.querySelectorAll("nav a.item")];
const obs = new IntersectionObserver(es => es.forEach(e => {
  if (e.isIntersecting) links.forEach(a => a.classList.toggle("active", a.getAttribute("href") === "#" + e.target.id));
}), { rootMargin: "-30% 0px -60% 0px" });
document.querySelectorAll("section").forEach(s => obs.observe(s));
links[0].classList.add("active");
"""
COPY_CSS = ".btn .ok { display: none; } .btn.done .ok { display: inline; } .btn.done .go { display: none; }"

page = f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title_zh)}</title>
<style>{CSS}{COPY_CSS}</style>
</head>
<body>
<nav>
  <div class="brand">{L(title_zh, title_en)}</div>
  <div class="brand-sub">{esc(p.pipeline)}</div>
  <div class="seg"><button data-lang="zh">中文</button><button data-lang="en">EN</button></div>
  <a class="item" href="#overview">{L("概览", "Overview")}</a>
  <a class="item" href="#workflow">{L("流程与参数", "Workflow & parameters")}</a>
  <a class="item" href="#summary">{L("QC 汇总表", "QC summary")}</a>
  <div class="ext"><a href="{esc(p.multiqc)}" target="_blank">{L("MultiQC 报告", "MultiQC report")} ↗</a></div>
</nav>
<main>
<header class="top">
  <h1>{L(title_zh, title_en)}</h1>
  <div class="project">{esc(PROJECT)}<span class="path">{esc(str(OUTDIR))}</span></div>
  <div class="meta">{esc(p.pipeline)} {esc(p.version)} · {L("生成于", "generated")} {datetime.now():%Y-%m-%d %H:%M} · {L(f"{n} 个样本", f"{n} samples")}</div>
</header>

<section id="overview">
  <h2>{L("概览", "Overview")}</h2>
  <div class="tiles">{tiles_html}</div>
  <div class="card">
    <h3>{L("需关注的样本", "Samples to review")}</h3>
    {attention}
  </div>
</section>

<section id="workflow">
  <h2>{L("流程与参数", "Workflow & parameters")}</h2>
  <div class="card"><div class="flow">{svg}</div></div>
  <div class="card">
    <h3>{L("本次参数", "Parameters of this run")}</h3>
    <table class="kv"><tbody>{settings_html}</tbody></table>
    <p class="muted">{L("链方向、文库类型、UMI 为各样本实际使用的值（样本表中填写的优先，留空时用 config 的默认值）。完整配置见 `logs/run_<时间>/config.yaml`。",
                        "Strandedness, protocol and UMI are the values actually used per sample (sample sheet first, config defaults otherwise). Full configuration: `logs/run_<time>/config.yaml`.")}</p>
  </div>
</section>

<section id="summary">
  <h2>{L("QC 汇总表", "QC summary")}</h2>
  <div class="card">
    <div class="bar">
      <div class="muted legend">{L("reads 按单条 read 计（PE 为 R1 + R2），rate 以进入该步骤的 reads 为分母（%，见表头括号）。点击表头排序；鼠标停在表头上显示含义与参考值，停在标色单元格上显示原因。",
                                    "Reads are single reads (R1 + R2 for PE); rates are % of the reads entering each step (see header). Click a header to sort; hover a header for its definition and reference, or a colored cell for the reason.")}
        <span class="sw" style="background:var(--orange-bg);border:1px solid var(--orange)"></span>{L("需注意", "review")}
        <span class="sw" style="background:var(--red-bg);border:1px solid var(--red)"></span>{L("低于公开标准", "below standard")}</div>
      <div class="tools">
        <button class="btn ghost" id="cols">{L("显示列", "Columns")} ▾</button>
        <div class="colpick" id="colpick"><span class="all">{L("全部显示", "Show all")}</span></div>
        <button class="btn" id="copy"><span class="go">{L("复制表格", "Copy table")}</span><span class="ok">{L("已复制", "Copied")} ✓</span></button>
      </div>
    </div>
    <div class="scroll">{table_html}</div>
  </div>
  <div class="card digest">
    <ul class="points">{conclusion}</ul>
    <p class="muted">{L("参考来源", "References")}：{sources}</p>
  </div>
</section>
</main>
<script>{JS}</script>
</body>
</html>
"""
Path(snakemake.output[0]).write_text(page)
print(f"{n} 个样本，{len(flagged)} 个需关注", file=sys.stderr)
