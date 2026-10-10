#!/usr/bin/env python3
"""32_fig3.py — Fig3（2×2）残余由病因决定（核心）。
 a 各数据集 免疫/ECM 残余 q | b immune−ECM 森林图 | c 无偏基因集(GSE199392) | d 解混淆 GSE241894
"""
import os, re, gzip
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import _figstyle as S

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT = f"{B}/data"; OUT = f"{B}/out"; IMG = f"{B}/out"
S.apply_style(); CAT = S.CB; POINT = 'k'; EMPH = S.CB[3]
FC = 1.0
def log2(x): return np.log2(x + 1)

def load_gmt(fn):
    s = {}
    for line in open(os.path.join(DAT, fn)):
        p = line.rstrip("\n").split("\t"); s[p[0]] = [g for g in p[2:] if g]
    return s
mh = load_gmt("mh.all.v2023.2.Mm.symbols.gmt"); m5 = load_gmt("m5.go.bp.v2023.2.Mm.symbols.gmt")
IMM = {'INFLAMM': mh['HALLMARK_INFLAMMATORY_RESPONSE'], 'TNFA': mh['HALLMARK_TNFA_SIGNALING_VIA_NFKB'],
       'IL6': mh['HALLMARK_IL6_JAK_STAT3_SIGNALING'], 'COMPL': mh['HALLMARK_COMPLEMENT'],
       'INNATE': m5['GOBP_INNATE_IMMUNE_RESPONSE'], 'MAC_ACT': m5['GOBP_MACROPHAGE_ACTIVATION'],
       'LEUK': m5['GOBP_LEUKOCYTE_MIGRATION']}
ECM = {'EMT': mh['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'], 'COL_FIB': m5['GOBP_COLLAGEN_FIBRIL_ORGANIZATION'],
       'ECM_ASM': m5['GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY'], 'COL_BIO': m5['GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS']}
gi = pd.read_csv(f"{DAT}/Mus_musculus.gene_info.gz", sep="\t", low_memory=False); gi.columns = [c.lstrip('#') for c in gi.columns]
s2e = {}
for s_, xr in zip(gi['Symbol'], gi['dbXrefs'].fillna('')):
    for m in re.findall(r'Ensembl:(ENSMUSG\d+)', xr): s2e.setdefault(s_, m)

def qset(X, ctrl, fib, rev, panels):
    c0 = X[ctrl].mean(1); fc = X[fib].mean(1) - c0
    Sg = fc[fc.abs() >= FC].index
    q = ((X[rev].mean(1) - c0).reindex(Sg).abs() / fc.reindex(Sg).abs()).clip(0, None)
    e2s = {e: next((k for k, v in s2e.items() if v == e), '') for e in Sg}
    out = {}
    for c, gs in panels.items():
        hit = [e for e in Sg if e2s.get(e, '') in set(gs)]
        out[c] = q.reindex(hit).dropna().values if len(hit) >= 5 else None
    return out

# ---- 载入各数据 ----
X1 = np.log2(pd.read_excel(f"{DAT}/GSE199392_processed_gene_expression_matrix.xlsx").set_index("#ID").astype(float) + 1)
g392 = ('ctrl', 'fib', 'rev')
D392 = dict(ctrl=['P0-%d'%i for i in range(1,6)], fib=['P6-%d'%i for i in [1,2,3,5,6]], rev=['R6-%d'%i for i in range(1,6)])
dg = pd.read_csv(f"{DAT}/GSE135462_FPKM_recovered.csv").set_index("ens"); Xg = log2(dg[[f"S{i}" for i in range(1,14)]].astype(float))
D462 = dict(ctrl=['S11','S12','S13'], fib=['S5','S6'], rev=['S9','S10'])
dh = pd.read_csv(f"{DAT}/GSE241894_counts.csv.gz", index_col=0); Xh = log2(dh.astype(float))
D4894 = dict(ctrl=[c for c in Xh.columns if c.replace(' ','').startswith('CNTRL.')],
             fib=[c for c in Xh.columns if c.replace(' ','').startswith('CCl4.')],
             rev=[c for c in Xh.columns if 'CCl4+siHAS2' in c.replace(' ','')],
             rev4mu=[c for c in Xh.columns if 'CCl4+4MU' in c.replace(' ','')])
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
e2 = pd.read_csv(f"{DAT}/GSE196908_Exp2_rpkm.txt.gz", sep="\t", index_col=0).astype(float)
e2.columns = [id2g.get(str(c), 'NA') for c in e2.columns]; e2 = log2(e2)
D196 = dict(ctrl=['Chow vehicle w12'], fib=['DIO-NASH vehicle w12'], rev=['DIO-NASH Chow Reversal w12'])

SETS = [("GSE199392", X1, D392, "toxin"), ("GSE135462", Xg, D462, "toxin"),
        ("GSE241894", Xh, D4894, "toxin"), ("GSE196908", e2, D196, "metabolic")]
rows = []
for name, X, D, et in SETS:
    qi = qset(X, D['ctrl'], D['fib'], D['rev'], IMM)
    qe = qset(X, D['ctrl'], D['fib'], D['rev'], ECM)
    iv = np.concatenate([v for v in qi.values() if v is not None])
    ev = np.concatenate([v for v in qe.values() if v is not None])
    p = mannwhitneyu(iv, ev, alternative='two-sided').pvalue
    rows.append((name, et, iv.mean(), ev.mean(), p, iv, ev))
    print(f"{name:11s} {et:9s} immune={iv.mean():.3f} ECM={ev.mean():.3f} diff={iv.mean()-ev.mean():+.3f} p={p:.2e}")

def L(ax, l): S.plab(ax, l)
fig = plt.figure(figsize=(7.2, 4.6)); gs = fig.add_gridspec(2, 2, hspace=0.72, wspace=0.36,
                                                           left=0.09, right=0.975, top=0.94, bottom=0.16)
CIMM, CECM = CAT[1], CAT[0]   # immune=orange / ECM=blue

rng = np.random.default_rng(0)
def boot_ci(a, b, n=2000):
    d = np.array([rng.choice(a, a.size).mean() - rng.choice(b, b.size).mean() for _ in range(n)])
    return np.percentile(d, 2.5), np.percentile(d, 97.5)
info = []
for name, et, mi, me, _, iv, ev in rows:
    p2 = mannwhitneyu(iv, ev, alternative='two-sided').pvalue
    lo, hi = boot_ci(iv, ev)
    info.append((name, et, mi - me, lo, hi, p2))

# (a) 各数据集 immune/ECM 残余分布（小提琴）
ax = fig.add_subplot(gs[0, 0])
for i, r in enumerate(rows):
    for arr, off, c in [(r[5], -0.19, CIMM), (r[6], 0.19, CECM)]:
        vp = ax.violinplot([arr], positions=[i + off], widths=0.34, showextrema=False)
        for b in vp['bodies']:
            b.set_facecolor(c); b.set_edgecolor(c); b.set_alpha(0.6); b.set_linewidth(0.5)
        ax.scatter([i + off], [arr.mean()], s=15, c="white", edgecolor="k", lw=0.6, zorder=5)
ax.set_xticks(range(len(rows)))
ax.set_xticklabels(['199392\ntoxin', '135462\ntoxin', '241894\ntoxin', '196908\nmetabolic'], fontsize=5.8)
ax.set_ylabel("residual $q$"); ax.set_ylim(-0.05, 3.15)
ax.legend(handles=[plt.Line2D([], [], marker='o', ls='', ms=4, mfc=CIMM, mec='none', label='immune'),
                   plt.Line2D([], [], marker='o', ls='', ms=4, mfc=CECM, mec='none', label='ECM')],
          loc="upper right", fontsize=5.6, handletextpad=0.3)
L(ax, "a")

# (b) immune−ECM 差值森林（bootstrap 95%CI，双侧 p）
ax = fig.add_subplot(gs[0, 1])
for i, (name, et, diff, lo, hi, p2) in enumerate(info):
    c = CIMM if diff > 0 else CECM
    ax.errorbar([diff], [i], xerr=[[diff - lo], [hi - diff]], fmt='o', ms=4.5,
                color=c, ecolor=c, elinewidth=1, capsize=2, zorder=3)
    lab = "p<0.001" if p2 < 0.001 else f"p={p2:.3f}"
    ax.text(0.345, i, lab, va='center', ha='right', fontsize=5.8)
ax.axvline(0, c='k', lw=0.7)
ax.set_yticks(range(len(info))); ax.set_yticklabels([r[0] for r in info], fontsize=6)
ax.set_ylim(len(info) - 0.5, -0.5); ax.set_xlim(-0.22, 0.36)
ax.set_xlabel("immune $q$ $-$ ECM $q$")
L(ax, "b")

# (c) GSE199392 逐基因集（棒棒糖）
ax = fig.add_subplot(gs[1, 0])
qi = qset(X1, D392['ctrl'], D392['fib'], D392['rev'], IMM)
qe = qset(X1, D392['ctrl'], D392['fib'], D392['rev'], ECM)
labs = list(IMM) + list(ECM)
vals = [qi[k].mean() for k in IMM] + [qe[k].mean() for k in ECM]
cols = [CIMM]*len(IMM) + [CECM]*len(ECM)
for j, (v, c) in enumerate(zip(vals, cols)):
    ax.plot([0, v], [j, j], '-', c=c, lw=0.9, zorder=1)
    ax.scatter([v], [j], s=20, c=c, zorder=3)
    ax.text(v + 0.015, j, f"{v:.2f}", va='center', fontsize=5.4)
ax.axhline(len(IMM) - 0.5, c='0.6', lw=0.8, ls='--')
ax.set_yticks(range(len(vals))); ax.set_yticklabels(labs, fontsize=5.2)
for t, c in zip(ax.get_yticklabels(), cols): t.set_color(c)
ax.set_ylim(len(vals) - 0.5, -0.5); ax.set_xlim(0, 1.05)
ax.set_xlabel("residual $q$ (GSE199392)")
L(ax, "c")

# (d) GSE241894 三臂（斜率图：随逆转力度推进的残余）
ax = fig.add_subplot(gs[1, 1])
def qall(X, ctrl, fib, rev):
    c0 = X[ctrl].mean(1); fc = X[fib].mean(1) - c0
    Sg = fc[fc.abs() >= FC].index
    return ((X[rev].mean(1) - c0).reindex(Sg).abs() / fc.reindex(Sg).abs()).clip(0, None).values
def ci_mean(x, n=2000):
    d = np.array([rng.choice(x, x.size).mean() for _ in range(n)])
    return np.percentile(d, [2.5, 97.5])
arms = [("CCl4 only", D4894['fib']), ("+4MU", D4894['rev4mu']), ("+siHAS2", D4894['rev'])]
xi = np.arange(3); IM = []; EC = []; AL = []; IMci = []; ECci = []; ps = []
for nm, rev in arms:
    qi_ = qset(Xh, D4894['ctrl'], D4894['fib'], rev, IMM); qe_ = qset(Xh, D4894['ctrl'], D4894['fib'], rev, ECM)
    a_ = np.concatenate([v for v in qi_.values() if v is not None])
    b_ = np.concatenate([v for v in qe_.values() if v is not None])
    IM.append(a_.mean()); EC.append(b_.mean()); AL.append(np.median(qall(Xh, D4894['ctrl'], D4894['fib'], rev)))
    IMci.append(ci_mean(a_)); ECci.append(ci_mean(b_))
    ps.append(mannwhitneyu(a_, b_, alternative='two-sided').pvalue)
IM = np.array(IM); EC = np.array(EC); AL = np.array(AL); IMci = np.array(IMci).T; ECci = np.array(ECci).T
ax.plot(xi, AL, '--', c='0.55', lw=1.0, marker='^', ms=3.4, label='all signature', zorder=1)
ax.errorbar(xi, IM, yerr=[IM-IMci[0], IMci[1]-IM], fmt='o-', c=CIMM, ms=4, lw=1.1, capsize=2, label='immune', zorder=3)
ax.errorbar(xi, EC, yerr=[EC-ECci[0], ECci[1]-EC], fmt='s-', c=CECM, ms=4, lw=1.1, capsize=2, label='ECM', zorder=3)
for i, p2 in enumerate(ps):
    lab = "n.s." if p2 >= 0.05 else ("p<0.001" if p2 < 0.001 else f"p={p2:.3f}")
    ax.text(xi[i], max(IM[i], EC[i]) + 0.13, lab, ha='center', fontsize=5.6,
            color=(S.CB[3] if p2 < 0.05 else '0.4'))
ax.set_xticks(xi); ax.set_xticklabels([a[0] for a in arms], fontsize=6)
ax.set_ylim(0, 1.45); ax.set_xlim(-0.35, 2.5)
ax.set_xlabel("CCl4 model   (reversal efficacy $\\to$)")
ax.set_ylabel("residual $q$ (CCl4 model)")
ax.legend(loc="upper right", fontsize=5.6, handletextpad=0.4, labelspacing=0.3)
L(ax, "d")

S.save(fig, IMG, "Fig3")
print("Fig3 done")
