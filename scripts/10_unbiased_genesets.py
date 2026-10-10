#!/usr/bin/env python3
"""10_unbiased_genesets.py — 用手选基因集换成无偏库(MSigDB 小鼠 Hallmark + GO BP)复测
"免疫/巨噬轴撤不掉、ECM 撤得掉"的结论, 以排除 selection 偏倚.

指标(每数据集, 造模签名 S={|FC_fib|>=1} 内):
  late q   = q(末点)              (残余程度, 0=完全回落)
  rebound  = q(末点) - min(q(前期)) (触底后回升幅度)
零分布: 从 S 内按 |FC_fib| 五分位匹配抽同 size 的子集, 1000x, 单边 p.
对照集(应接近背景): OXPHOS / 脂肪酸代谢 / 细胞周期.
"""
import os, re, gzip
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT  = os.path.join(ROOT, "data")
rng = np.random.default_rng(20261010)
FC = 1.0

def load_gmt(fn):
    sets = {}
    for line in open(os.path.join(DAT, fn)):
        p = line.rstrip("\n").split("\t")
        sets[p[0]] = [g for g in p[2:] if g]
    return sets
m1 = load_gmt("mh.all.v2023.2.Mm.symbols.gmt")
m5 = load_gmt("m5.go.bp.v2023.2.Mm.symbols.gmt")

PANEL = {
    # 免疫/巨噬/炎症 (预期: retreat 慢 / 反弹强)
    'HALLMARK_INFLAMMATORY_RESPONSE': m1.get('HALLMARK_INFLAMMATORY_RESPONSE', []),
    'HALLMARK_TNFA_SIGNALING_VIA_NFKB': m1.get('HALLMARK_TNFA_SIGNALING_VIA_NFKB', []),
    'HALLMARK_IL6_JAK_STAT3_SIGNALING': m1.get('HALLMARK_IL6_JAK_STAT3_SIGNALING', []),
    'HALLMARK_COMPLEMENT': m1.get('HALLMARK_COMPLEMENT', []),
    'GOBP_INNATE_IMMUNE_RESPONSE': m5.get('GOBP_INNATE_IMMUNE_RESPONSE', []),
    'GOBP_INFLAMMATORY_RESPONSE': m5.get('GOBP_INFLAMMATORY_RESPONSE', []),
    'GOBP_MACROPHAGE_ACTIVATION': m5.get('GOBP_MACROPHAGE_ACTIVATION', []),
    'GOBP_LEUKOCYTE_MIGRATION': m5.get('GOBP_LEUKOCYTE_MIGRATION', []),
    # ECM / 胶原 (预期: 撤得掉)
    'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION': m1.get('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', []),
    'GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY': m5.get('GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY', []),
    'GOBP_COLLAGEN_FIBRIL_ORGANIZATION': m5.get('GOBP_COLLAGEN_FIBRIL_ORGANIZATION', []),
    'GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS': m5.get('GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS', []),
    # 对照 (预期: 背景)
    'HALLMARK_OXIDATIVE_PHOSPHORYLATION': m1.get('HALLMARK_OXIDATIVE_PHOSPHORYLATION', []),
    'HALLMARK_FATTY_ACID_METABOLISM': m1.get('HALLMARK_FATTY_ACID_METABOLISM', []),
    'GOBP_CELL_CYCLE': m5.get('GOBP_CELL_CYCLE', []),
}
for k, v in PANEL.items():
    assert len(v) > 0, f"empty set {k}"
print("面板集大小:", {k: len(v) for k, v in PANEL.items()})

gi = pd.read_csv(os.path.join(DAT, "Mus_musculus.gene_info.gz"), sep="\t", low_memory=False)
gi.columns = [c.lstrip('#') for c in gi.columns]
sym2ens = {}
for s, xr in zip(gi['Symbol'], gi['dbXrefs'].fillna('')):
    for m in re.findall(r'Ensembl:(ENSMUSG\d+)', xr):
        sym2ens.setdefault(s, m)

def qtable(X, ctrl, fib, times, last, G):
    c0 = X[G[ctrl]].mean(1)
    fc = X[G[fib]].mean(1) - c0
    S = fc[fc.abs() >= FC].index
    q = pd.DataFrame({t: ((X[G[t]].mean(1) - c0).reindex(S).abs() / fc.reindex(S).abs()).clip(0, None) for t in times})
    q['fc'] = fc.reindex(S).abs()
    early = [t for t in times if t != last]
    q['late_q'] = q[last]
    q['rebound'] = q[last] - q[early].min(1)
    return q

def test_set(q, genes, nperm=1000):
    ens = [sym2ens[g] for g in genes if g in sym2ens]
    hit = [e for e in ens if e in q.index]
    if len(hit) < 5:
        return None
    sub = q.loc[hit]
    bins = pd.qcut(q['fc'].rank(method='first'), 5, labels=False)
    need = bins.loc[hit].value_counts()
    out = {}
    for metric in ['late_q', 'rebound']:
        obs = sub[metric].mean()
        null = np.empty(nperm)
        for i in range(nperm):
            pick = []
            for b, k in need.items():
                pool = bins.index[bins == b]
                pick += list(rng.choice(pool, size=min(k, len(pool)), replace=False))
            null[i] = q.loc[pick, metric].mean()
        out[metric] = (obs, (null >= obs).mean(), null.mean(), null.std())
    out['n'] = len(hit)
    return out

def run(name, q):
    print(f"\n{'='*78}\n### {name}   签名内 |FC|>=1 n={len(q)}")
    rows = []
    for gs, genes in PANEL.items():
        r = test_set(q, genes)
        if r is None:
            rows.append((gs, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)); continue
        rows.append((gs, r['n'], r['late_q'][0], r['late_q'][1], r['rebound'][0], r['rebound'][1],
                     r['rebound'][3], r['rebound'][0] / max(1e-9, r['rebound'][2])))
    df = pd.DataFrame(rows, columns=['geneset', 'n', 'late_q', 'p_late', 'rebound', 'p_reb', 'reb_sd', 'reb_z'])
    df = df.sort_values('rebound', ascending=False)
    for _, x in df.iterrows():
        print(f"  {x['geneset'][:44]:44s} n={int(x['n']):4d}  late_q={x['late_q']:.3f}(p={x['p_late']:.3f})  "
              f"rebound={x['rebound']:+.3f}(p={x['p_reb']:.3f})")

# GSE199392
X1 = np.log2(pd.read_excel(os.path.join(DAT, "GSE199392_processed_gene_expression_matrix.xlsx")).set_index("#ID").astype(float) + 1)
G1 = {'ctrl': ['P0-%d' % i for i in range(1, 6)], 'fib': ['P6-%d' % i for i in [1, 2, 3, 5, 6]],
      'R1': ['R1-%d' % i for i in range(1, 6)], 'R3': ['R3-%d' % i for i in range(1, 6)], 'R6': ['R6-%d' % i for i in range(1, 6)]}
run("GSE199392 (CCl4, R1/R3/R6)", qtable(X1, 'ctrl', 'fib', ['R1', 'R3', 'R6'], 'R6', G1))

# GSE261310
d = pd.read_excel(os.path.join(DAT, "GSE261310_Matrix.GRCm38.xlsx")).rename(columns=lambda c: str(c))
d = d.rename(columns={d.columns[0]: 'sym', d.columns[1]: 'ens'})
d = d[['ens'] + [c for c in d.columns if re.match(r'^(WT|T8_\d|T8F1_\d|T8F2_\d|T8F4_\d)', c)]].set_index('ens')
X2 = np.log2(d.astype(float) + 1)
G2 = {'ctrl': [c for c in X2.columns if c.startswith('WT_')], 'fib': [c for c in X2.columns if re.match(r'T8_\d', c)],
      'F1': [c for c in X2.columns if re.match(r'T8F1_\d', c)], 'F2': [c for c in X2.columns if re.match(r'T8F2_\d', c)],
      'F4': [c for c in X2.columns if re.match(r'T8F4_\d', c)]}
run("GSE261310 (TAA, F1/F2/F4)", qtable(X2, 'ctrl', 'fib', ['F1', 'F2', 'F4'], 'F4', G2))

# GSE135462
dg = pd.read_csv(os.path.join(DAT, "GSE135462_FPKM_recovered.csv")).set_index("ens")
X3 = np.log2(dg[[f"S{i}" for i in range(1, 14)]].astype(float) + 1)
G3 = {'ctrl': ['S11', 'S12', 'S13'], 'fib': ['S5', 'S6'], 'R1W': ['S7', 'S8'], 'R3W': ['S9', 'S10']}
run("GSE135462 (CCl4, R1W/R3W)", qtable(X3, 'ctrl', 'fib', ['R1W', 'R3W'], 'R3W', G3))
print("\nDONE")
