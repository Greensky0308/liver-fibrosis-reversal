#!/usr/bin/env python3
"""08_third_dataset.py — 控制项③: 用第三套独立数据检验"巨噬模块残余/反弹 > ECM".

注: GEO 无第三套 >=3 恢复时点的肝纤维化时程数据(已检索). 退而用"峰值 vs 恢复"
数据集(有正常对照)检验模块特异残余 q=|FC_恢复|/|FC_峰值| (0=完全回落,1=未回落).

纳入:
  GSE173961  TAA 8w -> Recovery, 对照 NS (n=3/组)
  GSE135462  CCl4 -> R1W/R3W 两个恢复点 (若可用)
"""
import os, re
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT  = os.path.join(ROOT, "data"); CACHE = os.path.join(ROOT, "out")
rng = np.random.default_rng(20261010)
FC = 1.0

gi = pd.read_csv(os.path.join(DAT, "Mus_musculus.gene_info.gz"), sep="\t", low_memory=False)
gi.columns = [c.lstrip('#') for c in gi.columns]
sym = {}
for s, xr in zip(gi['Symbol'], gi['dbXrefs'].fillna('')):
    for m in re.findall(r'Ensembl:(ENSMUSG\d+)', xr):
        sym[m] = s
sym2ens = {}
for e, s in sym.items():
    sym2ens.setdefault(s, e)

MOD = {
    'ECM_structural': ['Col1a1','Col1a2','Col3a1','Col4a1','Col4a2','Col5a1','Col5a2','Col6a1','Col6a2','Col6a3','Col8a1','Eln','Fn1','Lama4','Lamb1','Lamb2','Lamc1','Vcan','Fmod','Bgn','Dcn','Sparc','Postn'],
    'MMP_TIMP':       ['Mmp2','Mmp7','Mmp9','Mmp12','Mmp13','Mmp14','Mmp19','Timp1','Timp2','Timp3'],
    'crosslink_LOX':  ['Lox','Loxl1','Loxl2','Loxl3','Loxl4'],
    'HSC_activation': ['Acta2','Des','Pdgfrb','Tgfb1','Tgfbr1','Spp1','Tagln','Ctgf','Ccn2','Timp1'],
    'macrophage':     ['Cd68','Cd163','Mrc1','Lyve1','Ccl2','Ccl6','Ccl7','Adgre1','Trem2','Cd36','Csf1r','Msr1','Spp1'],
    'lipid_metab':    ['Cd36','Fgf21','Cidea','Cidec','Mogat1','Gpat3','Plin4','G0s2'],
}

def module_test(name, X, ctrl, fib, recs):
    """X: log2 expr (genes x samples); recs: {label: [cols]}"""
    print("\n" + "=" * 72)
    print(f"### {name}  对照组n={len(ctrl)} 纤维化n={len(fib)}")
    c0 = X[ctrl].mean(1)
    fc_fib = X[fib].mean(1) - c0
    sig = fc_fib[fc_fib.abs() >= FC].index
    print(f"纤维化签名 |FC|>=1: {len(sig)}")
    ens2sym = {e: sym.get(e, '') for e in sig}
    modens = {m: [e for e in sig if ens2sym.get(e, '') in set(gs)] for m, gs in MOD.items()}
    annotated = set(sum(modens.values(), []))
    for label, cols in recs.items():
        q = (X[cols].mean(1) - c0).reindex(sig).abs() / fc_fib.reindex(sig).abs()
        q = q.clip(0, None)
        q = pd.DataFrame({'q': q, 'fc': fc_fib.reindex(sig).abs()})
        rest = q[~q.index.isin(annotated)]
        print(f"\n  [{label}] 签名整体 q 中位={q['q'].median():.3f}")
        print(f"    {'module':16s} {'n':>3s} {'mean_q':>8s} {'p_vs_rest':>10s} {'null_mean':>10s} {'p_null':>8s}")
        bins = pd.qcut(q['fc'].rank(method='first'), 5, labels=False)
        for m in MOD:
            hit = modens[m]
            if len(hit) < 3:
                continue
            sub = q.loc[hit]
            try:
                _, p = mannwhitneyu(sub['q'], rest['q'], alternative='greater')
            except ValueError:
                p = np.nan
            need = bins.loc[hit].value_counts()
            null = []
            for _ in range(1000):
                pick = []
                for b, k in need.items():
                    pool = bins.index[bins == b]
                    pick += list(rng.choice(pool, size=min(k, len(pool)), replace=False))
                null.append(q.loc[pick, 'q'].mean())
            null = np.array(null)
            print(f"    {m:16s} {len(hit):3d} {sub['q'].mean():8.3f} {p:10.2e} "
                  f"{null.mean():10.3f} {(null >= sub['q'].mean()).mean():8.3f}")

# ---- GSE173961 ----
d = pd.read_csv(os.path.join(DAT, "GSE173961_mRNA_FPKM.txt.gz"), sep="\t").set_index("gene_id")
d = d[[c for c in d.columns if re.match(r'^(NS|Re|TAA)\d+$', c)]].astype(float)
X = np.log2(d + 1)
module_test("GSE173961 (TAA 8w -> Recovery)",
            X, ['NS1','NS2','NS3'], ['TAA1','TAA2','TAA3'], {'Recovery': ['Re1','Re2','Re3']})
print("\nDONE")
