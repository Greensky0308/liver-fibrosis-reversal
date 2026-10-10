#!/usr/bin/env python3
"""12_human_GSE106737.py — 人源锚: NASH 患者 RYGB 减重前后配对活检 (41 对).
问: 人类减重性消退中, ECM/胶原 vs 免疫/炎症 哪一轴撤得慢? (小鼠饮食撤除预测: ECM 滞后)
参照对比: 同数据集 patients_* 组 NASH(n=16) vs NONASH(n=13) 定义"纤维化方向 FC_ref".
配对变化: delta = followup - baseline (41 对).
消退比例 f = -delta/FC_ref (对 NASH 上调基因; 1=完全回到 NONASH 水平, 0=没动).
"""
import os, re, pickle
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT  = os.path.join(ROOT, "data"); CACHE = os.path.join(ROOT, "out")
FC = 1.0

# ---- 基因集 (人类) ----
def load_gmt(fn):
    s = {}
    for line in open(os.path.join(DAT, fn)):
        p = line.rstrip("\n").split("\t"); s[p[0]] = [g for g in p[2:] if g]
    return s
hs = load_gmt("h.all.v2023.2.Hs.symbols.gmt"); go = load_gmt("c5.go.bp.v2023.2.Hs.symbols.gmt")
IMMUNE = {'HALLMARK_INFLAMMATORY_RESPONSE': hs['HALLMARK_INFLAMMATORY_RESPONSE'],
          'HALLMARK_TNFA_NFKB': hs['HALLMARK_TNFA_SIGNALING_VIA_NFKB'],
          'HALLMARK_IL6_JAK_STAT3': hs['HALLMARK_IL6_JAK_STAT3_SIGNALING'],
          'HALLMARK_COMPLEMENT': hs['HALLMARK_COMPLEMENT'],
          'GOBP_INNATE_IMMUNE': go['GOBP_INNATE_IMMUNE_RESPONSE'],
          'GOBP_INFLAMMATORY': go['GOBP_INFLAMMATORY_RESPONSE'],
          'GOBP_MACROPHAGE_ACTIVATION': go['GOBP_MACROPHAGE_ACTIVATION'],
          'GOBP_LEUKOCYTE_MIGRATION': go['GOBP_LEUKOCYTE_MIGRATION']}
ECM = {'HALLMARK_EMT': hs['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'],
       'GOBP_COLLAGEN_FIBRIL_ORG': go['GOBP_COLLAGEN_FIBRIL_ORGANIZATION'],
       'GOBP_ECM_ASSEMBLY': go['GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY'],
       'GOBP_COLLAGEN_BIOSYNTH': go['GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS']}

# ---- 样本结构 ----
import gzip
for line in gzip.open(os.path.join(DAT, "GSE106737_series_matrix.txt.gz"), 'rt', errors='ignore'):
    if line.startswith('!Sample_title'):
        titles = [x.strip('"') for x in line.rstrip('\n').split('\t')[1:]]; break
hdr, rows = pickle.load(open(os.path.join(DAT, "gse106737_raw.pkl"), "rb"))
assert len(hdr) - 1 == len(titles), (len(hdr), len(titles))

# ---- probe -> symbol ----
m = pd.read_csv(os.path.join(DAT, "gpl16686_map.tsv"), sep="\t")
m = m.dropna(subset=["SYMBOL"]).drop_duplicates("PROBEID")
p2s = dict(zip(m["PROBEID"].astype(str), m["SYMBOL"]))
Xp = pd.DataFrame({hdr[0]: hdr[0]} and {r[0]: r[1:] for r in rows}).T.astype(float)
Xp.index.name = "probe"
Xp["sym"] = [p2s.get(str(i), "") for i in Xp.index]
Xg = Xp[Xp["sym"] != ""].groupby("sym").mean()          # 多探针取均值
Xg.columns = titles
print("gene-level matrix:", Xg.shape, " range", round(float(Xg.values.min()),2), round(float(Xg.values.max()),2))

# ---- 配对 ----
pairs = {}
for i, t in enumerate(titles):
    mm = re.match(r'Pair-(\d+)-', t)
    if mm:
        pairs.setdefault(mm.group(1), {})['baseline' if t.endswith('baseline') else 'followup'] = t
complete = {k: v for k, v in pairs.items() if len(v) == 2}
print("完整配对:", len(complete))
B = [v['baseline'] for v in complete.values()]; F = [v['followup'] for v in complete.values()]

# ---- 配对模块变化 (队列内, 无外部参照) ----
# 每基因配对 t 统计量 t_g = mean(delta)/se(delta); 负=改善后下调(纤维化基因回落), 接近0=动得少
Bn = Xg[B].values; Fn = Xg[F].values
D = Fn - Bn                                    # n=41 x genes
tstat = D.mean(1) / (D.std(1, ddof=1) / np.sqrt(D.shape[1]) + 1e-9)
tstat = pd.Series(tstat, index=Xg.index)
mdelta = pd.Series(D.mean(1), index=Xg.index)
print(f"\n配对 t: 中位={tstat.median():.3f}  |t|>2 基因数={(tstat.abs()>2).sum()}")

def modstat(panel, tag):
    print(f"\n[{tag}] per-set: mean t (负=整体下调/回落), mean Δ(log2)")
    out = {}
    for name, genes in panel.items():
        hit = [g for g in genes if g in tstat.index]
        if len(hit) < 5:
            print(f"   {name:28s} n={len(hit)} (跳过)"); continue
        out[name] = (len(hit), tstat.reindex(hit).mean(), mdelta.reindex(hit).mean())
        print(f"   {name:28s} n={len(hit):4d}  mean_t={tstat.reindex(hit).mean():+.3f}  meanΔ={mdelta.reindex(hit).mean():+.3f}")
    return out

mi = modstat(IMMUNE, 'IMMUNE'); me = modstat(ECM, 'ECM')

def pool(panel):
    gs = [g for genes in panel.values() for g in genes if g in tstat.index]
    return tstat.reindex(sorted(set(gs)))
ti, te = pool(IMMUNE), pool(ECM)
print(f"\n合并: 免疫 n={len(ti)} mean_t={ti.mean():+.3f} | ECM n={len(te)} mean_t={te.mean():+.3f}")
print(f"  MWU ECM>免疫(ECM下调更少=撤得慢) p={mannwhitneyu(te, ti, alternative='greater').pvalue:.3e}")
print(f"  MWU 免疫>ECM(免疫撤得慢)         p={mannwhitneyu(ti, te, alternative='greater').pvalue:.3e}")
print(f"  Δ: 免疫 mean={mdelta.reindex(sorted(set([g for gg in IMMUNE.values() for g in gg if g in mdelta.index]))).mean():+.3f}"
      f"  ECM mean={mdelta.reindex(sorted(set([g for gg in ECM.values() for g in gg if g in mdelta.index]))).mean():+.3f}")
print("\nDONE")
