#!/usr/bin/env python3
"""31_fig2.py — Fig2（2×2）逆转非单调。house style + 图型多样。
 a 折线(失调度随时间) | b 散点(触底 vs 末点 q, 每基因) | c 柱状(动力学分类) | d 热图(代表基因 q)
"""
import os, re
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import _figstyle as S

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT = f"{B}/data"; OUT = f"{B}/out"; IMG = f"{B}/out"
S.apply_style(); CB = S.CB
FC = 1.0
C1, C2 = CB[0], CB[1]

X1 = np.log2(pd.read_excel(f"{DAT}/GSE199392_processed_gene_expression_matrix.xlsx").set_index("#ID").astype(float) + 1)
G1 = {'P0': ['P0-%d'%i for i in range(1,6)], 'P6': ['P6-%d'%i for i in [1,2,3,5,6]],
      'R1': ['R1-%d'%i for i in range(1,6)], 'R3': ['R3-%d'%i for i in range(1,6)], 'R6': ['R6-%d'%i for i in range(1,6)]}
cnt = lambda df, ctrl, cols: int(((df[cols].mean(1) - df[ctrl].mean(1)).abs() >= FC).sum())
n1 = [cnt(X1, G1['P0'], G1[t]) for t in ['P6','R1','R3','R6']]
d = pd.read_excel(f"{DAT}/GSE261310_Matrix.GRCm38.xlsx"); d = d.rename(columns={d.columns[0]:'sym', d.columns[1]:'ens'})
d = d[['ens'] + [c for c in d.columns if re.match(r'^(WT|T8_\d|T8F1_\d|T8F2_\d|T8F4_\d)', str(c))]].set_index('ens')
X2 = np.log2(d.astype(float) + 1)
G2 = {'WT':[c for c in X2.columns if c.startswith('WT_')], 'T8':[c for c in X2.columns if re.match(r'T8_\d',c)],
      'F1':[c for c in X2.columns if re.match(r'T8F1_\d',c)], 'F2':[c for c in X2.columns if re.match(r'T8F2_\d',c)],
      'F4':[c for c in X2.columns if re.match(r'T8F4_\d',c)]}
n2 = [cnt(X2, G2['WT'], G2[t]) for t in ['T8','F1','F2','F4']]
kin = pd.read_csv(f"{OUT}/kinetics.csv", index_col=0)
kin = kin[kin['symbol'].astype(bool)]
print("n1", n1, "n2", n2, "genes", len(kin))

fig = plt.figure(figsize=(7.2, 4.6))
gs = fig.add_gridspec(2, 2, hspace=0.66, wspace=0.34, left=0.09, right=0.975, top=0.94, bottom=0.15)

# a 折线: 失调度随时间
ax = fig.add_subplot(gs[0, 0])
ax.plot([0,1,3,6], n1, "-o", ms=3.2, lw=1.1, c=C1, label="CCl4")
ax.plot([0,1,2,4], n2, "-s", ms=3.2, lw=1.1, c=C2, label="TAA")
ax.scatter([3],[min(n1)], s=55, facecolor="none", edgecolor=CB[3], zorder=5, linewidths=1.1)
ax.scatter([2],[min(n2)], s=55, facecolor="none", edgecolor=CB[3], zorder=5, linewidths=1.1)
ax.set_xlabel("weeks after withdrawal"); ax.set_ylabel("genes |log$_2$ FC| $\\geq$ 1")
ax.set_ylim(0, max(n1+n2)*1.22); ax.legend(loc="upper right", fontsize=5.4, handlelength=1.4, labelspacing=0.3)
ax.text(2.5, 3400, "trough", fontsize=5.8, color='0.35', ha='center', va='center')
ax.annotate("", xy=(2.98, min(n1)), xytext=(2.83, 3280),
            arrowprops=dict(arrowstyle="->", color='0.6', lw=0.6))
ax.annotate("", xy=(2.02, min(n2)), xytext=(2.30, 3620),
            arrowprops=dict(arrowstyle="->", color='0.6', lw=0.6))
S.plab(ax, "a")

# b 散点: 触底(R3) vs 末点(R6) q
ax = fig.add_subplot(gs[0, 1])
q3 = kin['R3'].values; q6 = kin['R6'].values
reb = q6 > q3
ax.scatter(q3[~reb], q6[~reb], s=2.5, c='0.75', lw=0, rasterized=True, label='stable')
ax.scatter(q3[reb], q6[reb], s=2.5, c=CB[1], lw=0, rasterized=True, label='rebound')
ax.plot([0, 3.8], [0, 3.8], '--', c='k', lw=0.6)
ax.set_xlim(0, 3.8); ax.set_ylim(0, 3.8)
ax.set_xlabel("resolution $q$ at trough (3w)"); ax.set_ylabel("$q$ at end (6w)")
frac = 100 * reb.mean()
ax.text(0.03, 0.95, f"{frac:.0f}% of genes rebound", transform=ax.transAxes,
        fontsize=5.8, va='top', bbox=dict(facecolor='white', edgecolor='none', pad=0.6))
ax.legend(loc="upper right", markerscale=3, fontsize=5.4, handletextpad=0.4, labelspacing=0.3)
S.plab(ax, "b")

# c 折线: 各类别 q 中位轨迹（动力学形状）
ax = fig.add_subplot(gs[1, 0])
order = ['fast','late','slow','rebound','persistent']
cols = [CB[5], CB[2], CB[4], CB[1], CB[0]]
for oc, cl in zip(cols, order):
    sub = kin[kin['class'] == cl][['R1','R3','R6']]
    med = sub.median()
    ax.plot([1,3,6], med.values, "-o", ms=2.6, lw=1.1, c=oc, label=f"{cl} ({len(sub)})")
ax.set_xlim(0.6, 6.4); ax.set_ylim(0, 1.32)
ax.set_xlabel("weeks after withdrawal"); ax.set_ylabel("median resolution $q$")
ax.legend(loc="upper right", fontsize=5.6, ncol=2, handlelength=1.2, labelspacing=0.2, columnspacing=0.9)
S.plab(ax, "c")

# d 哑铃: 代表基因 q 触底(3w)→末点(6w)
ax = fig.add_subplot(gs[1, 1])
k = kin.set_index('symbol')
grp = {g: CB[0] for g in ['Col1a1','Col3a1','Postn','Timp1']}
grp.update({g: CB[1] for g in ['Cd36','Cd68','Msr1','Ccl7']})
genes = [g for g in grp if g in k.index]
ys = np.arange(len(genes))[::-1]
for y, g in zip(ys, genes):
    q3, q6 = k.loc[g, 'R3'], k.loc[g, 'R6']
    c = grp[g]
    ax.plot([q3, q6], [y, y], "-", c=c, lw=1.0, alpha=0.75, zorder=1)
    ax.scatter([q3], [y], s=16, facecolor="white", edgecolor=c, lw=0.9, zorder=2)
    ax.scatter([q6], [y], s=17, c=c, zorder=3)
ax.axhline(3.5, c='0.85', lw=0.6, zorder=0)
ax.text(0.72, 5.5, "ECM / structural", color=CB[0], fontsize=5.4, ha='right', va='center')
ax.text(1.52, 1.5, "immune / macrophage", color=CB[1], fontsize=5.4, ha='right', va='center')
ax.set_yticks(ys); ax.set_yticklabels(genes, fontsize=6)
ax.set_xlim(0, 1.55); ax.set_ylim(-0.6, len(genes)-0.4)
ax.set_xlabel("resolution $q$")
h = [plt.Line2D([],[], marker='o', ls='', ms=3.6, mfc='white', mec='0.35', mew=0.9, label='3w trough'),
     plt.Line2D([],[], marker='o', ls='', ms=3.6, mfc='0.35', mec='0.35', label='6w end')]
ax.legend(handles=h, loc="upper right", fontsize=5.4, handletextpad=0.4, labelspacing=0.3)
S.plab(ax, "d")

S.save(fig, IMG, "Fig2")
print("Fig2 done")
