#!/usr/bin/env python3
"""09_GSE135462.py — 控制项③ 第三套: CCl4 进展(1/3/6W) + 逆转(R1W,R3W).
FPKM 从 GSE135462_m2.xlsx 抢救重建(sharedStrings 尾部损坏, 容错解压).
检验: (1) 逆转是否非单调 (q(R3W) vs q(R1W)); (2) 巨噬模块残余/反弹 vs ECM (匹配零分布).
"""
import os, re
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT  = os.path.join(ROOT, "data")
rng = np.random.default_rng(20261010)
FC = 1.0

gi = pd.read_csv(os.path.join(DAT, "Mus_musculus.gene_info.gz"), sep="\t", low_memory=False)
gi.columns = [c.lstrip('#') for c in gi.columns]
sym = {}
for s, xr in zip(gi['Symbol'], gi['dbXrefs'].fillna('')):
    for m in re.findall(r'Ensembl:(ENSMUSG\d+)', xr):
        sym[m] = s

MOD = {
    'ECM_structural': ['Col1a1','Col1a2','Col3a1','Col4a1','Col4a2','Col5a1','Col5a2','Col6a1','Col6a2','Col6a3','Col8a1','Eln','Fn1','Lama4','Lamb1','Lamb2','Lamc1','Vcan','Fmod','Bgn','Dcn','Sparc','Postn'],
    'MMP_TIMP':       ['Mmp2','Mmp7','Mmp9','Mmp12','Mmp13','Mmp14','Mmp19','Timp1','Timp2','Timp3'],
    'crosslink_LOX':  ['Lox','Loxl1','Loxl2','Loxl3','Loxl4'],
    'HSC_activation': ['Acta2','Des','Pdgfrb','Tgfb1','Tgfbr1','Spp1','Tagln','Ctgf','Ccn2','Timp1'],
    'macrophage':     ['Cd68','Cd163','Mrc1','Lyve1','Ccl2','Ccl6','Ccl7','Adgre1','Trem2','Cd36','Csf1r','Msr1','Spp1'],
    'lipid_metab':    ['Cd36','Fgf21','Cidea','Cidec','Mogat1','Gpat3','Plin4','G0s2'],
}

d = pd.read_csv(os.path.join(DAT, "GSE135462_FPKM_recovered.csv")).set_index("ens")
fpkm = d[[f"S{i}" for i in range(1, 14)]].astype(float)
X = np.log2(fpkm + 1)
G = {'ctrl': ['S11', 'S12', 'S13'], 'fib': ['S5', 'S6'],
     'R1W': ['S7', 'S8'], 'R3W': ['S9', 'S10']}
c0 = X[G['ctrl']].mean(1)
fc_fib = X[G['fib']].mean(1) - c0
sig = fc_fib[fc_fib.abs() >= FC].index
print(f"GSE135462 CCl4: 签名 |FC_fib|>=1 = {len(sig)}")
qdf = {}
for t in ['R1W', 'R3W']:
    qdf[t] = ((X[G[t]].mean(1) - c0).reindex(sig).abs() / fc_fib.reindex(sig).abs()).clip(0, None)
qdf = pd.DataFrame(qdf)
print(f"签名内 q 中位: R1W={qdf['R1W'].median():.3f}  R3W={qdf['R3W'].median():.3f}  "
      f"-> {'非单调(再恶化)' if qdf['R3W'].median() > qdf['R1W'].median() else '单调下降'}")

ens2sym = {e: sym.get(e, '') for e in sig}
modens = {m: [e for e in sig if ens2sym.get(e, '') in set(gs)] for m, gs in MOD.items()}
annotated = set(sum(modens.values(), []))
rest = qdf[~qdf.index.isin(annotated)]
bins = pd.qcut(fc_fib.reindex(sig).abs().rank(method='first'), 5, labels=False)
for t in ['R1W', 'R3W']:
    print(f"\n[{t}] 模块残余 q")
    print(f"  {'module':16s} {'n':>3s} {'mean_q':>8s} {'p_vs_rest':>10s} {'null':>8s} {'p_null':>8s}")
    for m in MOD:
        hit = modens[m]
        if len(hit) < 3: continue
        sub = qdf.loc[hit, t]
        try: _, p = mannwhitneyu(sub, rest[t], alternative='greater')
        except ValueError: p = np.nan
        need = bins.loc[hit].value_counts()
        null = []
        for _ in range(1000):
            pick = []
            for b, k in need.items():
                pool = bins.index[bins == b]
                pick += list(rng.choice(pool, size=min(k, len(pool)), replace=False))
            null.append(qdf.loc[pick, t].mean())
        null = np.array(null)
        print(f"  {m:16s} {len(hit):3d} {sub.mean():8.3f} {p:10.2e} {null.mean():8.3f} "
              f"{(null >= sub.mean()).mean():8.3f}")
print("\nDONE")
