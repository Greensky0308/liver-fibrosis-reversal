#!/usr/bin/env python3
"""30_fig1.py — Fig1（2×2）三种逆转模式并排。house style + 图型多样。
 a 柱状(各数据集签名大小) | b 棒棒糖(逆转比例) | c 小提琴+散点(可逆比) | d UpSet 点阵(签名重叠)
"""
import os, re, gzip
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import _figstyle as S

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT = f"{B}/data"; OUT = f"{B}/out"; IMG = f"{B}/out"
S.apply_style(); CB = S.CB
FC = 1.0
def log2(x): return np.log2(x + 1)
def sm_id2group(path):
    desc = tit = None
    for line in gzip.open(path, 'rt', errors='ignore'):
        f = line.rstrip('\n').split('\t')
        if line.startswith('!Sample_description'): desc = f[1:]
        if line.startswith('!Sample_title'): tit = f[1:]
    desc = [d.strip('"') for d in desc]
    tit = [re.sub(r'^Exp\s*2\s*', '', re.sub(r'\s*Rep\s*\d+\s*$', '', t.strip('"'))) for t in tit]
    return dict(zip(desc, tit))
id2g = sm_id2group(f"{DAT}/GSE196908_series_matrix.txt.gz")
def load_exp(f):
    d = pd.read_csv(f"{DAT}/{f}", sep="\t", index_col=0).astype(float)
    d.columns = [id2g.get(str(c), 'NA') for c in d.columns]; return log2(d)
x1 = log2(pd.read_excel(f"{DAT}/GSE199392_processed_gene_expression_matrix.xlsx").set_index("#ID").astype(float))
e1 = load_exp("GSE196908_rpkm.txt.gz"); e2 = load_exp("GSE196908_Exp2_rpkm.txt.gz")
d = pd.read_excel(f"{DAT}/GSE261310_Matrix.GRCm38.xlsx"); d = d.rename(columns={d.columns[0]: 'sym', d.columns[1]: 'ens'})
d = d[['ens'] + [c for c in d.columns if re.match(r'^(WT|T8_\d|T8F1_\d|T8F2_\d|T8F4_\d)', str(c))]].set_index('ens')
X2 = log2(d.astype(float))
Xg = log2(pd.read_csv(f"{DAT}/GSE135462_FPKM_recovered.csv").set_index("ens")[[f"S{i}" for i in range(1,14)]].astype(float))
Xh = log2(pd.read_csv(f"{DAT}/GSE241894_counts.csv.gz", index_col=0).astype(float))
def sig(df, ctrl, fib):
    A = df[fib].values; Bv = df[ctrl].values
    fc = A.mean(1) - Bv.mean(1); p = ttest_ind(A, Bv, axis=1, equal_var=False).pvalue
    return set(np.asarray(df.index)[np.abs(fc) >= FC])
cA = ['P0-%d' % i for i in range(1, 6)]
DS = [("199392\ntoxin", x1, cA, ['P6-%d' % i for i in [1,2,3,5,6]], CB[0]),
      ("261310\ntoxin", X2, [c for c in X2.columns if c.startswith('WT_')], [c for c in X2.columns if re.match(r'T8_\d', c)], CB[0]),
      ("135462\ntoxin", Xg, ['S11','S12','S13'], ['S5','S6'], CB[0]),
      ("241894\ntoxin", Xh, [c for c in Xh.columns if c.replace(' ','').startswith('CNTRL.')], [c for c in Xh.columns if c.replace(' ','').startswith('CCl4.')], CB[0]),
      ("196908\ndiet", e2, ['Chow vehicle w12'], ['DIO-NASH vehicle w12'], CB[2]),
      ("196908\ndrug", e1, ['Chow+Vehicle 12w'], ['DIO-NASH Vehicle 12w'], CB[1])]
sizes = [len(sig(df, c, f)) for _, df, c, f, _ in DS]
colr = [cc for *_, cc in DS]
sigA = sig(x1, cA, ['P6-%d' % i for i in [1,2,3,5,6]])
sigB = sig(e2, ['Chow vehicle w12'], ['DIO-NASH vehicle w12'])
sigC = sig(e1, ['Chow+Vehicle 12w'], ['DIO-NASH Vehicle 12w'])
def frac(df, sg, ctrl, res):
    r = (df[res].mean(1) - df[ctrl].mean(1)).reindex(list(sg)); return (r.abs() < FC/2).mean() * 100
revA_t = [frac(x1, sigA, cA, ['R%d-%d' % (w, i) for i in range(1, 6)]) for w in (1, 3, 6)]
revA = np.mean(revA_t)
revB = frac(e2, sigB, ['Chow vehicle w12'], ['DIO-NASH Chow Reversal w12'])
revC = np.mean([frac(e1, sigC, ['Chow+Vehicle 12w'], t) for t in ['Lanifibranor 12w', 'Semaglutide 12w']])
rr = pd.read_csv(f"{OUT}/reversibility_ratio.csv", index_col=0)
print("sizes", sizes, "rev", round(revA), round(revB), round(revC))

fig = plt.figure(figsize=(7.2, 4.6))
gs = fig.add_gridspec(2, 2, hspace=0.72, wspace=0.34, left=0.09, right=0.975, top=0.94, bottom=0.14)

# a 柱状: 各数据集签名大小
ax = fig.add_subplot(gs[0, 0])
ax.bar(np.arange(len(sizes)), sizes, width=0.66, color=colr)
for i, v in enumerate(sizes): ax.text(i, v + max(sizes)*0.02, str(v), ha='center', fontsize=5.2)
ax.set_xticks(np.arange(len(sizes))); ax.set_xticklabels([n for n, *_ in DS]); S.setx(ax, 5.6, 45, 'right')
ax.set_ylabel("signature genes"); ax.set_ylim(0, max(sizes)*1.22)
_h = [plt.Rectangle((0, 0), 1, 1, color=CB[0], label='toxin'),
      plt.Rectangle((0, 0), 1, 1, color=CB[2], label='diet'),
      plt.Rectangle((0, 0), 1, 1, color=CB[1], label='drug')]
ax.legend(handles=_h, loc='upper right', fontsize=5.4, handlelength=0.9, handleheight=0.9,
          labelspacing=0.3, borderpad=0.3)
S.plab(ax, "a")

# b 折线: 逆转比例随撤除时间 (毒物有时程)
ax = fig.add_subplot(gs[0, 1])
modes = ["toxin", "diet", "drug"]; rev = [revA, revB, revC]; modecol = [CB[0], CB[2], CB[1]]
ax.plot([1, 3, 6], revA_t, "-o", ms=4, lw=1.4, c=CB[0], label="toxin")
ax.axhline(revB, c=CB[2], lw=1.2, ls="--"); ax.text(8.7, revB, f"diet {revB:.0f}%", va="center", ha="center", fontsize=5.6, color=CB[2],
        bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
ax.axhline(revC, c=CB[1], lw=1.2, ls=":"); ax.text(8.7, revC, f"drug {revC:.0f}%", va="center", ha="center", fontsize=5.6, color=CB[1],
        bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
for x_, y_ in zip([1, 3, 6], revA_t):
    ax.text(x_, y_ + 2.6, f"{y_:.0f}", ha="center", va="bottom", fontsize=5.4, color=CB[0])
ax.set_xlabel("weeks after withdrawal"); ax.set_ylabel("signature reverted (%)")
ax.set_xlim(0.4, 10.8); ax.set_ylim(0, max(revA_t) * 1.35)
ax.legend(loc="upper right")
S.plab(ax, "b")

# c 小提琴 + 个体散点: 可逆比
ax = fig.add_subplot(gs[1, 0])
data = [rr['A_toxinW_r'].dropna().values, rr['B_dietW_r'].dropna().values, rr['C_drug_r'].dropna().values]
parts = ax.violinplot(data, positions=[0, 1, 2], widths=0.75, showextrema=False, showmedians=False)
rng = np.random.default_rng(0)
for k, (pc, cc) in enumerate(zip(parts['bodies'], modecol)):
    pc.set_facecolor(cc); pc.set_alpha(0.55); pc.set_edgecolor('k'); pc.set_linewidth(0.4)
    v = data[k]
    ax.plot([k-0.22, k+0.22], [np.median(v)]*2, c='k', lw=0.9)
ax.set_xticks([0, 1, 2]); ax.set_xticklabels(modes); ax.set_ylim(0, 1.1)
ax.axhline(0.5, c='0.8', lw=0.5, ls='--')
ax.set_ylabel("reversibility ratio $r$")
S.plab(ax, "c")

# d UpSet 点阵: 三模式签名重叠
sub = gs[1, 1].subgridspec(2, 1, height_ratios=[2.0, 1.2], hspace=0.08)
sets = {'toxin': sigA, 'diet': sigB, 'drug': sigC}
names = ['toxin', 'diet', 'drug']
inter = [('toxin',), ('diet',), ('drug',), ('toxin','diet'), ('toxin','drug'), ('diet','drug'), ('toxin','diet','drug')]
cnt = []
for combo in inter:
    s = set.intersection(*[sets[c] for c in combo])
    for o in names:
        if o not in combo: s = s - sets[o]
    cnt.append(len(s))
axT = fig.add_subplot(sub[0])
axT.bar(np.arange(len(inter)), cnt, width=0.62, color=[CB[0] if 'toxin' in c else (CB[2] if 'diet' in c else CB[1]) for c in inter])
for i, v in enumerate(cnt): axT.text(i, v + max(cnt)*0.03, str(v), ha='center', fontsize=5.2)
axT.set_xticks([]); axT.set_ylim(0, max(cnt)*1.25); axT.set_ylabel("genes", fontsize=6.2)
S.plab(axT, "d")
axM = fig.add_subplot(sub[1])
for r, o in enumerate(names):
    for i, combo in enumerate(inter):
        inS = o in combo
        axM.plot(i, r, 'o', ms=3.6, c=(CB[[0,2,1][r]] if inS else '0.85'), zorder=3)
for i, combo in enumerate(inter):
    rs = [r for r, o in enumerate(names) if o in combo]
    if len(rs) > 1: axM.plot([i, i], [min(rs), max(rs)], '-', c='0.4', lw=0.9, zorder=2)
axM.set_yticks(range(3)); axM.set_yticklabels(names, fontsize=5.8); axM.invert_yaxis()
axM.set_xticks(np.arange(len(inter))); axM.set_ylim(2.6, -0.6); axM.set_xlim(-0.6, len(inter)-0.4)
S.despine(axM)

S.save(fig, IMG, "Fig1")
print("Fig1 done")
