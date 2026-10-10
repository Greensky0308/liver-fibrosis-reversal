#!/usr/bin/env python3
"""07_specificity.py — 控制项①: 晚期"再激活"(反弹)是纤维化程序特异, 还是全转录组漂移?

方法修正: q(t)=|FC_t|/|FC_fib| 在小分母(|FC_fib|~0)下不稳定 -> 全部限制在
造模签名 S={|FC_fib|>=1} 内计算 (分母稳定).

三个层次的检验:
  L1 签名轨迹: q(R1/R3/R6) 中位 -> 是否非单调(触底后回升)   [稳健, 无小分母]
  L2 模块特异: 在 S 内, ECM/HSC/巨噬/... 的反弹 dq vs S 内其余基因 (MWU)
  L3 匹配零分布: 从 S 内按 |FC_fib| 五分位分层, 抽同 size 随机子集, 1000x -> p
判据: 若模块 dq 显著 > S 内其余/匹配零分布 => 再激活为该模块特异.
另附: 全转录组 |FC|>1 基因数随撤除时间的轨迹 (背景参考).
"""
import os, re
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT  = os.path.join(ROOT, "data"); CACHE = os.path.join(ROOT, "out")
rng = np.random.default_rng(20261009)
FC = 1.0

# ---------- 符号表 ----------
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
sym2ens = {}
for e, s in sym.items():
    sym2ens.setdefault(s, e)

def load_199392():
    X = np.log2(pd.read_excel(os.path.join(DAT, "GSE199392_processed_gene_expression_matrix.xlsx"))
                .set_index("#ID").astype(float) + 1)
    G = {'ctrl': ['P0-%d' % i for i in range(1, 6)], 'fib': ['P6-%d' % i for i in [1, 2, 3, 5, 6]],
         'R1': ['R1-%d' % i for i in range(1, 6)], 'R3': ['R3-%d' % i for i in range(1, 6)],
         'R6': ['R6-%d' % i for i in range(1, 6)]}
    return X, G, ['R1', 'R3', 'R6'], 'R6'

def load_261310():
    d = pd.read_excel(os.path.join(DAT, "GSE261310_Matrix.GRCm38.xlsx"))
    d = d.rename(columns={d.columns[0]: 'sym', d.columns[1]: 'ens'})
    keep = ['ens'] + [c for c in d.columns if re.match(r'^(WT|T8|T12|T8F)', str(c))]
    d = d[keep].set_index('ens')
    X = np.log2(d.astype(float) + 1)
    G = {'ctrl': [c for c in X.columns if re.match(r'WT_', c)],
         'fib':  [c for c in X.columns if re.match(r'T8_\d', c)],
         'F1':   [c for c in X.columns if re.match(r'T8F1_\d', c)],
         'F2':   [c for c in X.columns if re.match(r'T8F2_\d', c)],
         'F4':   [c for c in X.columns if re.match(r'T8F4_\d', c)]}
    return X, G, ['F1', 'F2', 'F4'], 'F4'

def analyze(name, X, G, tps, last):
    print("\n" + "=" * 72)
    print(f"### {name}   对照n={len(G['ctrl'])} 造模n={len(G['fib'])} 撤除n={[len(G[t]) for t in tps]} 末点={last}")
    c0 = X[G['ctrl']].mean(1)
    fc_fib = (X[G['fib']].mean(1) - c0)
    fc_t = {t: (X[G[t]].mean(1) - c0) for t in tps}
    sig = fc_fib[fc_fib.abs() >= FC].index
    print(f"造模签名 |FC_fib|>=1: {len(sig)}")

    # ---- 全转录组 |FC|>1 基因数轨迹 (背景参考) ----
    print("\n[全转录组 |FC|>1 基因数]  " + "  ".join(
        f"{t}={int((fc_t[t].abs() > FC).sum())}" for t in tps))

    # ---- L1: 签名内 q 轨迹 (分母稳定) ----
    q = pd.DataFrame({t: (fc_t[t].reindex(sig).abs() / fc_fib.reindex(sig).abs()) for t in tps})
    qmed = {t: q[t].median() for t in tps}
    trough = min(tps, key=lambda t: qmed[t])
    print("[L1 签名内 q 中位] " + "  ".join(f"q({t})={qmed[t]:.3f}" for t in tps)
          + f"   -> 触底于 {trough}, 末点{'高于' if qmed[last] > qmed[trough] else '未高于'}触底 "
          f"({'非单调(再激活)' if qmed[last] > qmed[trough] else '单调'})")
    early = [t for t in tps if t != last]
    q['qmin'] = q[early].min(1)
    q['dq']   = q[last] - q['qmin']
    q['fc']   = fc_fib.reindex(sig).abs()
    print(f"        S 内整体 dq: mean={q['dq'].mean():+.3f} median={q['dq'].median():+.3f} "
          f"frac(dq>0.3)={np.mean(q['dq'] > 0.3):.1%}")
    rho, prho = spearmanr(q['fc'], q['dq'])
    print(f"        Spearman(|FC_fib|, dq) rho={rho:+.3f} p={prho:.1e}")

    # ---- L2: 模块 dq vs S 内其余 ----
    ens2sym = {e: sym.get(e, '') for e in sig}
    modens = {}
    for m, gs in MOD.items():
        modens[m] = [e for e in sig if ens2sym.get(e, '') in set(gs)]
    annotated = set(sum(modens.values(), []))
    rest = q[~q.index.isin(annotated)]
    print("\n[L2 模块反弹 dq vs S内其余基因]")
    print(f"    {'module':16s} {'n':>3s} {'mean_dq':>8s} {'frac>0.3':>8s} {'mean_qmin':>10s} {'p_vs_rest':>10s}")
    for m in MOD:
        hit = modens[m]
        if not hit:
            continue
        s = q.loc[hit]
        try:
            _, p = mannwhitneyu(s['dq'], rest['dq'], alternative='greater')
        except ValueError:
            p = np.nan
        print(f"    {m:16s} {len(hit):3d} {s['dq'].mean():+8.3f} {np.mean(s['dq']>0.3):8.1%} "
              f"{s['qmin'].mean():10.2f} {p:10.2e}")

    # ---- L3: |FC_fib| 分层匹配零分布 ----
    print("\n[L3 匹配零分布 (S 内按 |FC_fib| 五分位分层抽同 size)]")
    bins = pd.qcut(q['fc'].rank(method='first'), 5, labels=False)
    for m in MOD:
        hit = modens[m]
        if len(hit) < 3:
            continue
        obs = q.loc[hit, 'dq'].mean()
        # 每层需要抽的个数 = 模块内落在该层的个数
        need = bins.loc[hit].value_counts()
        null = []
        for _ in range(1000):
            pick = []
            for b, k in need.items():
                pool = bins.index[bins == b]
                pick += list(rng.choice(pool, size=min(k, len(pool)), replace=False))
            null.append(q.loc[pick, 'dq'].mean())
        null = np.array(null)
        print(f"    {m:16s} obs={obs:+.3f}  null mean={null.mean():+.3f} sd={null.std():.3f}  "
              f"p={(null >= obs).mean():.3f}")
    return q

q1 = analyze("GSE199392 (CCl4, 撤除 1/3/6 周)", *load_199392())
q2 = analyze("GSE261310 (TAA, 撤除 1/2/4 周)", *load_261310())
q1.assign(symbol=[sym.get(e, '') for e in q1.index]).to_csv(os.path.join(CACHE, "specificity_GSE199392.csv"))
q2.assign(symbol=[sym.get(e, '') for e in q2.index]).to_csv(os.path.join(CACHE, "specificity_GSE261310.csv"))
print("\nDONE")
