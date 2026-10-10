#!/usr/bin/env python3
"""15_mcp_counter.py — 用 MCP-counter 正规签名(人源→小鼠大小写)估细胞丰度.
对称检验: 毒物模型(CCl4)的"免疫残余"是否=免疫细胞残留;
          代谢模型(GAN-diet)的"ECM残余"是否=成纤维细胞残留.
细胞群: T / CD8T / Cytotoxic / NK / B / Monocytic / MyeloidDC / Neutrophils / Endothelial / Fibroblasts.
"""
import os, re, gzip
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT  = os.path.join(ROOT, "data")
FC = 1.0

# ---- MCP-counter 签名 (人 -> 鼠大小写) ----
sig = {}
for line in open(os.path.join(DAT, "mcpc.txt")):
    p = [x.strip('"') for x in line.rstrip("\n").split("\t")]
    if len(p) < 2 or p[0] == "HUGO symbols": continue
    sig.setdefault(p[1], []).append(p[0].capitalize())
print("MCP-counter 群:", {k: len(v) for k, v in sig.items()})

gi = pd.read_csv(os.path.join(DAT, "Mus_musculus.gene_info.gz"), sep="\t", low_memory=False)
gi.columns = [c.lstrip('#') for c in gi.columns]
sym2ens = {}
for s, xr in zip(gi['Symbol'], gi['dbXrefs'].fillna('')):
    for m in re.findall(r'Ensembl:(ENSMUSG\d+)', xr): sym2ens.setdefault(s, m)

def score(X, name):
    print(f"\n{'='*72}\n### {name}  (样本 {X.shape[1]})")
    rows = {}
    for pop, genes in sig.items():
        hit = [sym2ens[g] for g in genes if g in sym2ens and sym2ens[g] in X.index]
        rows[pop] = X.loc[hit].mean(0) if hit else pd.Series(0.0, index=X.columns)
    sc = pd.DataFrame(rows)
    return sc

def traj(sc, G, order):
    out = pd.DataFrame({t: sc.loc[G[t]].mean() for t in order})
    print(out.round(2).to_string())

# ---- GSE199392 (CCl4) ----
X1 = np.log2(pd.read_excel(os.path.join(DAT, "GSE199392_processed_gene_expression_matrix.xlsx")).set_index("#ID").astype(float) + 1)
G1 = {'ctrl': ['P0-%d'%i for i in range(1,6)], 'fib': ['P6-%d'%i for i in [1,2,3,5,6]],
      'R1': ['R1-%d'%i for i in range(1,6)], 'R3': ['R3-%d'%i for i in range(1,6)], 'R6': ['R6-%d'%i for i in range(1,6)]}
s1 = score(X1, "GSE199392 CCl4")
print("丰度轨迹 (ctrl/fib/R1/R3/R6):"); traj(s1, G1, ['ctrl','fib','R1','R3','R6'])

# ---- GSE196908 (GAN-diet) ----
def sm_id2group(path):
    desc = tit = None
    for line in gzip.open(path, 'rt', errors='ignore'):
        f = line.rstrip('\n').split('\t')
        if line.startswith('!Sample_description'): desc = f[1:]
        if line.startswith('!Sample_title'): tit = f[1:]
    desc = [d.strip('"') for d in desc]
    tit = [re.sub(r'^Exp\s*2\s*', '', re.sub(r'\s*Rep\s*\d+\s*$', '', t.strip('"'))) for t in tit]
    return dict(zip(desc, tit))
id2g = sm_id2group(os.path.join(DAT, "GSE196908_series_matrix.txt.gz"))
def load_exp(f):
    d = pd.read_csv(os.path.join(DAT, f), sep="\t", index_col=0).astype(float)
    d.columns = [id2g.get(str(c), 'NA') for c in d.columns]; return np.log2(d + 1)
e1 = load_exp("GSE196908_rpkm.txt.gz"); e2 = load_exp("GSE196908_Exp2_rpkm.txt.gz")
G_e2 = {'ctrl': ['Chow vehicle w12'], 'fib': ['DIO-NASH vehicle w12'], 'rev': ['DIO-NASH Chow Reversal w12']}
s2 = score(e2, "GSE196908 Exp2 GAN-diet (饮食撤除)")
print("丰度轨迹 (Chow/NASH-veh/Chow-Reversal):"); traj(s2, G_e2, ['ctrl','fib','rev'])
G_e1 = {'ctrl': ['Chow+Vehicle 12w'], 'fib': ['DIO-NASH Vehicle 12w'], 'lani': ['Lanifibranor 12w'], 'sema': ['Semaglutide 12w']}
s1b = score(e1, "GSE196908 Exp1 GAN-diet (药物)")
print("丰度轨迹 (Chow/NASH-veh/Lani/Sema):"); traj(s1b, G_e1, ['ctrl','fib','lani','sema'])

# ---- 关键: 残余 q 与细胞丰度残差的关系 ----
def residual_vs_abundance(X, G, times, last, sc, pop):
    c0 = X[G['ctrl']].mean(1); fc = X[G['fib']].mean(1) - c0
    S = fc[fc.abs() >= FC].index
    q = ((X[G[last]].mean(1) - c0).reindex(S).abs() / fc.reindex(S).abs()).clip(0, None)
    # 丰度: ctrl->fib->last
    ab = {t: sc[pop].loc[G[t]].mean() for t in ['ctrl','fib',last]}
    print(f"  [{pop}] 丰度 ctrl={ab['ctrl']:.2f} fib={ab['fib']:.2f} {last}={ab[last]:.2f}  "
          f"-> 残余比例=(last-ctrl)/(fib-ctrl)={ (ab[last]-ab['ctrl'])/(ab['fib']-ab['ctrl']+1e-9):.2f}")
    return q

print(f"\n{'='*72}\n### 细胞丰度残余 (0=完全回落,1=没回落) — 毒物 vs 代谢")
print("GSE199392 CCl4:")
for pop in ['Monocytic','T cells','NK cells','B lineage','Neutrophils','Fibroblasts','Endothelial cells']:
    if pop in s1.columns: residual_vs_abundance(X1, G1, ['R1','R3','R6'], 'R6', s1, pop)
print("GSE196908 Exp2 GAN-diet:")
for pop in ['Monocytic','T cells','NK cells','B lineage','Neutrophils','Fibroblasts','Endothelial cells']:
    if pop in s2.columns: residual_vs_abundance(e2, G_e2, ['rev'], 'rev', s2, pop)
print("\nDONE")
