#!/usr/bin/env python3
"""33_fig4.py — Fig4（2×2）人源锚: NASH RYGB 前后 41 对配对活检。
 a 设计示意(文本) | b 各基因集 mean t | c 逐基因配对 t 分布 | d 人 vs 小鼠 方向一致
"""
import os, re, gzip, pickle
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import _figstyle as S

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT = f"{B}/data"; OUT = f"{B}/out"; IMG = f"{B}/out"
S.apply_style(); CAT = S.CB; POINT = 'k'; EMPH = S.CB[3]
CI, CE = CAT[0], CAT[2]     # immune / ECM

def load_gmt(fn):
    s = {}
    for line in open(os.path.join(DAT, fn)):
        p = line.rstrip("\n").split("\t"); s[p[0]] = [g for g in p[2:] if g]
    return s
hs = load_gmt("h.all.v2023.2.Hs.symbols.gmt"); go = load_gmt("c5.go.bp.v2023.2.Hs.symbols.gmt")
IMM = {'INFLAMM': hs['HALLMARK_INFLAMMATORY_RESPONSE'], 'TNFA': hs['HALLMARK_TNFA_SIGNALING_VIA_NFKB'],
       'IL6': hs['HALLMARK_IL6_JAK_STAT3_SIGNALING'], 'COMPL': hs['HALLMARK_COMPLEMENT'],
       'INNATE': go['GOBP_INNATE_IMMUNE_RESPONSE'], 'MAC_ACT': go['GOBP_MACROPHAGE_ACTIVATION'],
       'LEUK': go['GOBP_LEUKOCYTE_MIGRATION']}
ECM = {'EMT': hs['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'], 'COL_FIB': go['GOBP_COLLAGEN_FIBRIL_ORGANIZATION'],
       'ECM_ASM': go['GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY'], 'COL_BIO': go['GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS']}

# ---- 载入人源 ----
for line in gzip.open(f"{DAT}/GSE106737_series_matrix.txt.gz", 'rt', errors='ignore'):
    if line.startswith('!Sample_title'):
        titles = [x.strip('"') for x in line.rstrip('\n').split('\t')[1:]]; break
hdr, rows = pickle.load(open(os.path.join(DAT, "gse106737_raw.pkl"), "rb"))
m = pd.read_csv(os.path.join(DAT, "gpl16686_map.tsv"), sep="\t").dropna(subset=["SYMBOL"]).drop_duplicates("PROBEID")
p2s = dict(zip(m["PROBEID"].astype(str), m["SYMBOL"]))
Xp = pd.DataFrame({r[0]: r[1:] for r in rows}).T.astype(float)
Xp["sym"] = [p2s.get(str(i), "") for i in Xp.index]
Xg = Xp[Xp["sym"] != ""].groupby("sym").mean(); Xg.columns = titles
pairs = {}
for t in titles:
    mm = re.match(r'Pair-(\d+)-', t)
    if mm: pairs.setdefault(mm.group(1), {})['baseline' if t.endswith('baseline') else 'followup'] = t
comp = {k: v for k, v in pairs.items() if len(v) == 2}
Bc = [v['baseline'] for v in comp.values()]; F = [v['followup'] for v in comp.values()]
D = Xg[F].values - Xg[Bc].values
tstat = pd.Series(D.mean(1) / (D.std(1, ddof=1) / np.sqrt(D.shape[1]) + 1e-9), index=Xg.index)

def setv(panels):
    out = {}
    for c, gs in panels.items():
        hit = [g for g in gs if g in tstat.index]
        if len(hit) >= 5: out[c] = tstat.reindex(hit).values
    return out
vi = setv(IMM); ve = setv(ECM)
alli = np.concatenate(list(vi.values())); alle = np.concatenate(list(ve.values()))
p_all = mannwhitneyu(alle, alli, alternative='greater').pvalue
print(f"human: immune mean_t={alli.mean():.3f} (n={len(alli)})  ECM mean_t={alle.mean():.3f} (n={len(alle)})  p(ECM>immune)={p_all:.2e}")

fig = plt.figure(figsize=(7.2, 4.6)); gs = fig.add_gridspec(2, 2, hspace=0.72, wspace=0.34, left=0.09, right=0.975, top=0.94, bottom=0.15)

# (a) 每患者 免疫/ECM 模块评分 前后配对
ax = fig.add_subplot(gs[0, 0])
imm_gs = [g for gs in IMM.values() for g in gs if g in Xg.index]
ecm_gs = [g for gs in ECM.values() for g in gs if g in Xg.index]
si = Xg.loc[imm_gs].mean(0); se = Xg.loc[ecm_gs].mean(0)
si = (si - si.mean()) / si.std(); se = (se - se.mean()) / se.std()
bi = [si[v['baseline']] for v in comp.values()]; fi = [si[v['followup']] for v in comp.values()]
be = [se[v['baseline']] for v in comp.values()]; fe = [se[v['followup']] for v in comp.values()]
for b_, f_ in zip(bi, fi): ax.plot([0, 1], [b_, f_], '-', color=CI, lw=0.4, alpha=0.35)
for b_, f_ in zip(be, fe): ax.plot([0, 1], [b_, f_], '-', color=CE, lw=0.4, alpha=0.35)
ax.plot([0, 1], [np.mean(bi), np.mean(fi)], '-o', color=CI, lw=1.6, ms=3.5, label='immune')
ax.plot([0, 1], [np.mean(be), np.mean(fe)], '-o', color=CE, lw=1.6, ms=3.5, label='ECM')
ax.set_xticks([0, 1]); ax.set_xticklabels(["before", "after\nRYGB"])
ax.set_ylabel("module score ($z$)"); ax.set_xlim(-0.25, 1.25)
ax.legend(loc="upper right", bbox_to_anchor=(1.06, 1.0), fontsize=5.4, handletextpad=0.4, labelspacing=0.3)
S.plab(ax, "a")

# (b) 各基因集 mean t
ax = fig.add_subplot(gs[0, 1])
labs = list(IMM) + list(ECM)
vals = [vi[k].mean() for k in IMM] + [ve[k].mean() for k in ECM]
cols = [CI]*len(IMM) + [CE]*len(ECM)
ax.barh(np.arange(len(vals)), vals, color=cols)
ax.set_yticks(np.arange(len(vals))); ax.set_yticklabels(labs, fontsize=6.2); ax.invert_yaxis()
ax.axvline(0, c=POINT, lw=0.8)
ax.set_xlabel("mean paired $t$  (negative = down)"); ax.set_xlim(-1.9, 0.25)
ax.axhline(len(IMM)-0.5, c="0.6", lw=0.8, ls="--")
S.plab(ax, "b")

# (c) 逐基因配对 t 分布
ax = fig.add_subplot(gs[1, 0])
parts = ax.violinplot([alli, alle], positions=[0, 1], widths=0.8, showextrema=False, showmedians=True)
for pc, cc in zip(parts['bodies'], [CI, CE]):
    pc.set_facecolor(cc); pc.set_alpha(0.85); pc.set_edgecolor(POINT); pc.set_linewidth(0.6)
parts['cmedians'].set_color(POINT); parts['cmedians'].set_linewidth(1.2)
ax.axhline(0, c="0.8", lw=0.8, ls="--")
ax.set_xticks([0, 1]); ax.set_xticklabels(["immune", "ECM"])
ax.set_ylabel("per-gene paired $t$"); ax.set_ylim(-6.5, 5.5)
ax.text(0.5, 4.9, f"p = {p_all:.1e}", ha="center", fontsize=6.2, color=EMPH)
S.plab(ax, "c")

# (d) 逐基因集: 小鼠 q vs 人 mean_t 一致性
ax = fig.add_subplot(gs[1, 1])
def load_gmt_m(fn):
    s = {}
    for line in open(os.path.join(DAT, fn)):
        p = line.rstrip("\n").split("\t"); s[p[0]] = [g for g in p[2:] if g]
    return s
mhm = load_gmt_m("mh.all.v2023.2.Mm.symbols.gmt"); m5m = load_gmt_m("m5.go.bp.v2023.2.Mm.symbols.gmt")
NAME2M = {'INFLAMM': ('mh','HALLMARK_INFLAMMATORY_RESPONSE'), 'TNFA': ('mh','HALLMARK_TNFA_SIGNALING_VIA_NFKB'),
          'IL6': ('mh','HALLMARK_IL6_JAK_STAT3_SIGNALING'), 'COMPL': ('mh','HALLMARK_COMPLEMENT'),
          'INNATE': ('m5','GOBP_INNATE_IMMUNE_RESPONSE'), 'MAC_ACT': ('m5','GOBP_MACROPHAGE_ACTIVATION'),
          'LEUK': ('m5','GOBP_LEUKOCYTE_MIGRATION'), 'EMT': ('mh','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'),
          'COL_FIB': ('m5','GOBP_COLLAGEN_FIBRIL_ORGANIZATION'), 'ECM_ASM': ('m5','GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY'),
          'COL_BIO': ('m5','GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS')}
gi = pd.read_csv(f"{DAT}/Mus_musculus.gene_info.gz", sep="\t", low_memory=False); gi.columns = [c.lstrip('#') for c in gi.columns]
s2e = {}
for s_, xr in zip(gi['Symbol'], gi['dbXrefs'].fillna('')):
    for mm_ in re.findall(r'Ensembl:(ENSMUSG\d+)', xr): s2e.setdefault(s_, mm_)
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
e2.columns = [id2g.get(str(c), 'NA') for c in e2.columns]; e2 = np.log2(e2 + 1)
c0 = e2[['Chow vehicle w12']].mean(1); fc = e2[['DIO-NASH vehicle w12']].mean(1) - c0
Sg = fc[fc.abs() >= 1].index
q = ((e2[['DIO-NASH Chow Reversal w12']].mean(1) - c0).reindex(Sg).abs() / fc.reindex(Sg).abs()).clip(0, None)
e2s = {e: next((k for k, v in s2e.items() if v == e), '') for e in Sg}
xm, ym, cm = [], [], []
for k in IMM.keys() | ECM.keys():
    src, gname = NAME2M[k]
    gs = (mhm if src == 'mh' else m5m)[gname]
    hit = [e for e in Sg if e2s.get(e, '') in set(gs)]
    if len(hit) < 5 or k not in (set(IMM) | set(ECM)): continue
    human_v = (vi.get(k, np.array([])).mean() if k in vi else ve.get(k, np.array([np.nan])).mean())
    xm.append(q.reindex(hit).dropna().mean()); ym.append(human_v)
    cm.append(CI if k in IMM else CE)
ax.scatter(xm, ym, s=70, c=cm, zorder=3)
ax.axhline(0, c="0.85", lw=0.7); ax.axvline(0.5, c="0.85", lw=0.7)
ax.set_xlabel("mouse residual $q$ (GSE196908)"); ax.set_ylabel("human mean paired $t$")
ax.set_xlim(0.25, 0.85); ax.set_ylim(-1.9, 0.4)
ax.legend(handles=[plt.Line2D([], [], marker='o', ls='', ms=4.5, mfc=CI, mec='none', label='immune'),
                   plt.Line2D([], [], marker='o', ls='', ms=4.5, mfc=CE, mec='none', label='ECM')],
          loc="upper right", fontsize=5.4, handletextpad=0.3, labelspacing=0.3)
S.plab(ax, "d")

S.save(fig, IMG, "Fig4")
print("Fig4 done")
