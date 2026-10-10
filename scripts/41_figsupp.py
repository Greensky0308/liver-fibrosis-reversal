#!/usr/bin/env python3
"""41_figsupp.py — 补充图 S1-S3。
S1 细胞去卷积：免疫残余 = 免疫细胞未退走（MCP-counter + marker 丰度校正）
S2 单细胞 GSE295293 三方法交叉验证（聚类 / ComBat / 逐细胞打分）
S3 反弹特异性（|FC| 分层匹配零分布）+ 第三数据集 GSE135462
"""
import os, re, gzip
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import _figstyle as S

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT = f"{B}/data"; OUT = f"{B}/out"; IMG = f"{B}/out"
S.apply_style(); CB = S.CB; POINT = 'k'
FC = 1.0
rng = np.random.default_rng(20261010)
CI, CE = CB[1], CB[0]          # immune=orange / ECM=blue
CTX, CMT = CB[0], CB[2]        # toxin / metabolic


def load_gmt(fn):
    s = {}
    for line in open(f"{DAT}/{fn}"):
        p = line.rstrip("\n").split("\t"); s[p[0]] = [g for g in p[2:] if g]
    return s
mh = load_gmt("mh.all.v2023.2.Mm.symbols.gmt"); m5 = load_gmt("m5.go.bp.v2023.2.Mm.symbols.gmt")
IMM_SET = {'INFLAMM': mh['HALLMARK_INFLAMMATORY_RESPONSE'], 'TNFA': mh['HALLMARK_TNFA_SIGNALING_VIA_NFKB'],
           'IL6': mh['HALLMARK_IL6_JAK_STAT3_SIGNALING'], 'COMPL': mh['HALLMARK_COMPLEMENT'],
           'INNATE': m5['GOBP_INNATE_IMMUNE_RESPONSE'], 'MAC_ACT': m5['GOBP_MACROPHAGE_ACTIVATION'],
           'LEUK': m5['GOBP_LEUKOCYTE_MIGRATION']}
ECM_SET = {'EMT': mh['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'], 'COL_FIB': m5['GOBP_COLLAGEN_FIBRIL_ORGANIZATION'],
           'ECM_ASM': m5['GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY'], 'COL_BIO': m5['GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS']}

gi = pd.read_csv(f"{DAT}/Mus_musculus.gene_info.gz", sep="\t", low_memory=False); gi.columns = [c.lstrip('#') for c in gi.columns]
sym = {}; sym2ens = {}
for s_, xr in zip(gi['Symbol'], gi['dbXrefs'].fillna('')):
    for m in re.findall(r'Ensembl:(ENSMUSG\d+)', xr):
        sym[m] = s_; sym2ens.setdefault(s_, m)

CELL = {'Mac_Kupffer': ['Cd68', 'Csf1r', 'Adgre1', 'Mrc1', 'Cd163', 'Lyz2', 'Msr1', 'Trem2', 'Vsig4', 'Marco', 'Timd4', 'Clec4f', 'Fcgr1'],
        'Monocyte': ['Ccr2', 'Ly6c2', 'Chil3', 'S100a8', 'S100a9', 'Fcnb', 'Ms4a6c'],
        'T_cell': ['Cd3e', 'Cd3d', 'Cd3g', 'Cd8a', 'Cd4', 'Trbc2', 'Lck', 'Thy1'],
        'NK': ['Nkg7', 'Klrb1c', 'Ncr1', 'Klrd1', 'Prf1'], 'B_cell': ['Cd79a', 'Cd79b', 'Ms4a1', 'Cd19', 'Ighm'],
        'Neutrophil': ['Ly6g', 'Mpo', 'Elane', 'Retnlg', 'S100a8', 'S100a9'],
        'HSC_mesench': ['Acta2', 'Des', 'Pdgfrb', 'Gfap', 'Gli1', 'Col1a1', 'Col3a1'],
        'Hepatocyte': ['Alb', 'Apoa1', 'Apoa2', 'Ttr', 'Serpina1a', 'Cyp2e1', 'Asgr1', 'Fga'],
        'Endothelial': ['Pecam1', 'Cdh5', 'Kdr', 'Stab2', 'Clec4g', 'Lyve1']}
IMMUNE_CELLS = ['Mac_Kupffer', 'Monocyte', 'T_cell', 'NK', 'B_cell', 'Neutrophil']
IDENTITY = set(sum((CELL[c] for c in IMMUNE_CELLS), []))


def sm_trunc(t):
    return re.sub(r'^Exp\s*2\s*', '', re.sub(r'\s*Rep\s*\d+\s*$', '', t.strip('"')))


def load_196908(f):
    desc = tit = None
    for line in gzip.open(f"{DAT}/GSE196908_series_matrix.txt.gz", 'rt', errors='ignore'):
        fx = line.rstrip('\n').split('\t')
        if line.startswith('!Sample_description'): desc = fx[1:]
        if line.startswith('!Sample_title'): tit = fx[1:]
    id2g = dict(zip([d.strip('"') for d in desc], [sm_trunc(t) for t in tit]))
    d = pd.read_csv(f"{DAT}/{f}", sep="\t", index_col=0).astype(float)
    d.columns = [id2g.get(str(c), 'NA') for c in d.columns]
    return np.log2(d + 1)


# ============================== Fig S1 ==============================
# MCP-counter 签名 (人 -> 鼠大小写)
mcpsig = {}
for line in open(os.path.join(DAT, "mcpc.txt")):
    p = [x.strip('"') for x in line.rstrip("\n").split("\t")]
    if len(p) < 2 or p[0] == "HUGO symbols": continue
    mcpsig.setdefault(p[1], []).append(p[0].capitalize())


def mcp_score(X):
    rows = {}
    for pop, genes in mcpsig.items():
        hit = [sym2ens[g] for g in genes if g in sym2ens and sym2ens[g] in X.index]
        rows[pop] = X.loc[hit].mean(0) if hit else pd.Series(0.0, index=X.columns)
    return pd.DataFrame(rows)


X1 = np.log2(pd.read_excel(f"{DAT}/GSE199392_processed_gene_expression_matrix.xlsx").set_index("#ID").astype(float) + 1)
G1 = {'ctrl': ['P0-%d' % i for i in range(1, 6)], 'fib': ['P6-%d' % i for i in [1, 2, 3, 5, 6]],
      'R1': ['R1-%d' % i for i in range(1, 6)], 'R3': ['R3-%d' % i for i in range(1, 6)], 'R6': ['R6-%d' % i for i in range(1, 6)]}
e2 = load_196908("GSE196908_Exp2_rpkm.txt.gz")
G2 = {'ctrl': ['Chow vehicle w12'], 'fib': ['DIO-NASH vehicle w12'], 'rev': ['DIO-NASH Chow Reversal w12']}
sc1 = mcp_score(X1); sc2 = mcp_score(e2)
POPS = ['T cells', 'NK cells', 'B lineage', 'Monocytic lineage', 'Fibroblasts', 'Endothelial cells']


def resid(sc, G, last):
    return [((sc[pop].loc[G[last]].mean() - sc[pop].loc[G['ctrl']].mean()) /
             (sc[pop].loc[G['fib']].mean() - sc[pop].loc[G['ctrl']].mean() + 1e-9)) for pop in POPS]


r_tox = resid(sc1, G1, 'R6'); r_met = resid(sc2, G2, 'rev')

# marker 去卷积: 免疫 vs ECM 残余 q (校正前/后), GSE199392
def deconv(name, X, G, times, last):
    c0 = X[G['ctrl']].mean(1)
    Xsym = {s_: X.loc[e] for s_, e in sym2ens.items() if e in X.index}
    sc = pd.DataFrame({ct: pd.concat([Xsym[g] for g in gs if g in Xsym], axis=1).mean(1)
                       for ct, gs in CELL.items() if any(g in Xsym for g in gs)})
    imm_sc = sc[IMMUNE_CELLS].sum(1)
    fc = X[G['fib']].mean(1) - c0; Sig = fc[fc.abs() >= FC].index
    ens2s = {e: sym.get(e, '') for e in Sig}
    imm_g = sorted({e for e in Sig if ens2s[e] in set(sum(IMM_SET.values(), []))})
    ecm_g = sorted({e for e in Sig if ens2s[e] in set(sum(ECM_SET.values(), []))})
    imm_id = [e for e in imm_g if ens2s[e] in IDENTITY]; imm_ef = [e for e in imm_g if ens2s[e] not in IDENTITY]
    q = ((X[G[last]].mean(1) - c0).reindex(Sig).abs() / fc.reindex(Sig).abs()).clip(0, None)
    allc = G['ctrl'] + G['fib'] + sum((G[t] for t in times), [])
    A = imm_sc.reindex(allc).values; A = (A - A.mean()) / (A.std() + 1e-9); A = np.column_stack([np.ones(len(A)), A])
    Xs = X.reindex(Sig)[allc].values
    beta, *_ = np.linalg.lstsq(A, Xs.T, rcond=None); resid = (Xs.T - A @ beta).T
    Xr = pd.DataFrame(resid, index=Sig, columns=allc); c0r = Xr[G['ctrl']].mean(1); fcr = Xr[G['fib']].mean(1) - c0r
    qr = ((Xr[G[last]].mean(1) - c0r).abs() / fcr.abs()).clip(0, None)
    return q, qr, imm_g, imm_ef, ecm_g


q1, qr1, ig1, ie1, eg1 = deconv("199392", X1, G1, ['R1', 'R3', 'R6'], 'R6')
cats = ['Immune\n(all)', 'Immune\nidentity', 'Immune\neffector', 'ECM']
before = [q1.reindex(ig1).mean(), q1.reindex([e for e in ig1 if sym.get(e, '') in IDENTITY]).mean() if False else q1.reindex(ig1).mean(), 0, q1.reindex(eg1).mean()]
# 明确重算 identity / effector
id_ens = [e for e in ig1 if sym.get(e, '') in IDENTITY]; ef_ens = [e for e in ig1 if sym.get(e, '') not in IDENTITY]
before = [q1.reindex(ig1).mean(), q1.reindex(id_ens).mean(), q1.reindex(ef_ens).mean(), q1.reindex(eg1).mean()]
after = [qr1.reindex(ig1).mean(), qr1.reindex(id_ens).mean(), qr1.reindex(ef_ens).mean(), qr1.reindex(eg1).mean()]
p_dec = mannwhitneyu(q1.reindex(ig1), q1.reindex(eg1), alternative='greater').pvalue
p_decr = mannwhitneyu(qr1.reindex(ig1), qr1.reindex(eg1), alternative='greater').pvalue
print("S1a MCP resid tox:", np.round(r_tox, 2), " met:", np.round(r_met, 2))
print("S1b before:", np.round(before, 3), "after:", np.round(after, 3), f"p={p_dec:.3f}->{p_decr:.3f}")

fig = plt.figure(figsize=(7.2, 3.0))
gs = fig.add_gridspec(1, 2, wspace=0.30, left=0.085, right=0.98, top=0.90, bottom=0.20)
ax = fig.add_subplot(gs[0, 0])
x = np.arange(len(POPS)); w = 0.38
ax.bar(x - w/2, r_tox, w, color=CTX, label='CCl4 (toxin)')
ax.bar(x + w/2, r_met, w, color=CMT, label='GAN-diet (metabolic)')
ax.axhline(1, c='0.5', lw=0.8, ls='--')   # 1 = no retreat
ax.set_xticks(x); ax.set_xticklabels(POPS, fontsize=5.6); S.setx(ax, 5.6, 40, 'right')
ax.set_ylabel("cell-abundance residual"); ax.set_ylim(0, max(r_tox + r_met) * 1.18 + 0.1)
ax.legend(loc='upper right', fontsize=5.4, handlelength=1.2, labelspacing=0.3)
ax.text(0.02, 0.96, "1 = no retreat", transform=ax.transAxes, fontsize=5.2, color='0.4', va='top')
S.plab(ax, "a")

ax = fig.add_subplot(gs[0, 1])
x = np.arange(len(cats)); w = 0.38
ax.bar(x - w/2, before, w, color='0.75', label='uncorrected')
ax.bar(x + w/2, after, w, color=CI, label='abundance-corrected')
ax.set_xticks(x); ax.set_xticklabels(cats, fontsize=5.4)
ax.set_ylabel("residual $q$ (GSE199392)"); ax.set_ylim(0, 1.05)
ax.legend(loc='upper right', fontsize=5.4, handlelength=1.2, labelspacing=0.3)
ax.text(0.02, 0.97, f"immune>ECM  p={p_dec:.3f} $\\to$ {p_decr:.3f}", transform=ax.transAxes,
        fontsize=5.2, color='0.3', va='top')
S.plab(ax, "b")
S.save(fig, IMG, "FigS1")


# ============================== Fig S2 ==============================
def prop_from(csv):
    a = pd.read_csv(f"{OUT}/{csv}")
    a = a[a['celltype'] != 'Unassigned']
    tot = a.groupby('group').size()
    tab = a.groupby(['group', 'celltype']).size().unstack(fill_value=0)
    return (tab.div(tot, axis=0) * 100)


methods = [("cluster", "GSE295293_celltype.csv"), ("ComBat", "GSE295293_celltype_combat.csv"),
           ("per-cell score", "GSE295293_celltype_percell.csv")]
CELLS = ['Kupffer/Mac', 'HSC/Mesench', 'Monocyte']; GROUPS = ['control', 'fibrosis', 'regression']
GCOL = {'control': CB[6], 'fibrosis': CB[1], 'regression': CB[0]}
PROPS = [prop_from(csv) for _, csv in methods]
GMAX = max(P[CELLS].max().max() for P in PROPS)

fig = plt.figure(figsize=(7.2, 2.8))
gs = fig.add_gridspec(1, 3, wspace=0.26, left=0.075, right=0.99, top=0.88, bottom=0.22)
for i, (mname, csv) in enumerate(methods):
    P = PROPS[i]
    ax = fig.add_subplot(gs[0, i])
    x = np.arange(len(CELLS)); w = 0.26
    for k, g in enumerate(GROUPS):
        v = [P.loc[g, c] if c in P.columns else 0 for c in CELLS]
        ax.bar(x + (k - 1) * w, v, w, color=GCOL[g], label=(g if i == 0 else None))
    ax.set_xticks(x); ax.set_xticklabels(CELLS, fontsize=5.4); S.setx(ax, 5.4, 30, 'right')
    ax.set_ylim(0, GMAX * 1.12)
    if i == 0:
        ax.set_ylabel("cells (%)"); ax.legend(loc='upper right', fontsize=5.0, handlelength=1.0, labelspacing=0.2)
    else:
        ax.set_yticklabels([])
    S.plab(ax, "abc"[i])
S.save(fig, IMG, "FigS2")


# ============================== Fig S3 ==============================
MOD = {'ECM_structural': ['Col1a1', 'Col1a2', 'Col3a1', 'Col4a1', 'Col4a2', 'Col5a1', 'Col5a2', 'Col6a1', 'Col6a2', 'Col6a3', 'Col8a1', 'Eln', 'Fn1', 'Lama4', 'Lamb1', 'Lamb2', 'Lamc1', 'Vcan', 'Fmod', 'Bgn', 'Dcn', 'Sparc', 'Postn'],
       'Macrophage': ['Cd68', 'Cd163', 'Mrc1', 'Lyve1', 'Ccl2', 'Ccl6', 'Ccl7', 'Adgre1', 'Trem2', 'Cd36', 'Csf1r', 'Msr1', 'Spp1']}


def module_null(csv, rng):
    q = pd.read_csv(f"{OUT}/{csv}", index_col=0)
    q = q[q['symbol'].astype(bool)]
    bins = pd.qcut(q['fc'].rank(method='first'), 5, labels=False)
    res = {}
    for m, gs in MOD.items():
        hit = q[q['symbol'].isin(set(gs))]
        obs = hit['dq'].mean()
        need = bins.loc[hit.index].value_counts()
        null = []
        for _ in range(1000):
            pick = []
            for b, k in need.items():
                pool = bins.index[bins == b]
                pick += list(rng.choice(pool, size=min(k, len(pool)), replace=False))
            null.append(q.loc[pick, 'dq'].mean())
        null = np.array(null)
        res[m] = (obs, null, (null >= obs).mean())
    return res


r1 = module_null("specificity_GSE199392.csv", rng)
r2 = module_null("specificity_GSE261310.csv", rng)

fig = plt.figure(figsize=(7.2, 3.0))
gs = fig.add_gridspec(1, 2, wspace=0.32, left=0.085, right=0.98, top=0.88, bottom=0.22)
ax = fig.add_subplot(gs[0, 0])
labels = ['GSE199392\nECM', 'GSE199392\nMacrophage', 'GSE261310\nECM', 'GSE261310\nMacrophage']
nulls = [r1['ECM_structural'][1], r1['Macrophage'][1], r2['ECM_structural'][1], r2['Macrophage'][1]]
obss = [r1['ECM_structural'][0], r1['Macrophage'][0], r2['ECM_structural'][0], r2['Macrophage'][0]]
ps = [r1['ECM_structural'][2], r1['Macrophage'][2], r2['ECM_structural'][2], r2['Macrophage'][2]]
cols = [CE, CI, CE, CI]
parts = ax.violinplot(nulls, positions=range(4), widths=0.8, showextrema=False, showmedians=False)
for pc, cc in zip(parts['bodies'], cols):
    pc.set_facecolor(cc); pc.set_alpha(0.28); pc.set_edgecolor(cc); pc.set_linewidth(0.7)
ax.axhline(0, c='0.7', lw=0.7, ls='--')
PY = 0.86
for i, (o, p, cc) in enumerate(zip(obss, ps, cols)):
    ax.scatter([i], [o], s=26, c=cc, zorder=5)
    ax.text(i, PY, f"p={p:.2f}", fontsize=5.2, color='0.25', va='bottom', ha='center')
ax.set_xticks(range(4)); ax.set_xticklabels(labels, fontsize=5.2)
ax.set_ylabel("module rebound $\\Delta q$"); ax.set_ylim(-0.35, 1.0)
ax.text(0.02, 0.97, "violin = matched null", transform=ax.transAxes, fontsize=5.2, color='0.4', va='top')
S.plab(ax, "a")

# GSE135462 第三数据集: 巨噬 vs ECM 残余 q
dg = pd.read_csv(f"{DAT}/GSE135462_FPKM_recovered.csv").set_index("ens")
Xg = np.log2(dg[[f"S{i}" for i in range(1, 14)]].astype(float) + 1)
Gg = {'ctrl': ['S11', 'S12', 'S13'], 'fib': ['S5', 'S6'], 'R1W': ['S7', 'S8'], 'R3W': ['S9', 'S10']}
c0 = Xg[Gg['ctrl']].mean(1); fcg = Xg[Gg['fib']].mean(1) - c0
Sg = fcg[fcg.abs() >= FC].index; es = {e: sym.get(e, '') for e in Sg}
mm = {'Macrophage': ['Cd68', 'Cd163', 'Mrc1', 'Lyve1', 'Ccl2', 'Ccl6', 'Ccl7', 'Adgre1', 'Trem2', 'Cd36', 'Csf1r', 'Msr1', 'Spp1'],
      'ECM': ['Col1a1', 'Col1a2', 'Col3a1', 'Col4a1', 'Col4a2', 'Col5a1', 'Col5a2', 'Col6a1', 'Col6a2', 'Col6a3', 'Col8a1', 'Eln', 'Fn1', 'Lama4', 'Lamb1', 'Lamb2', 'Lamc1', 'Vcan', 'Fmod', 'Bgn', 'Dcn', 'Sparc', 'Postn']}
qg = {t: ((Xg[Gg[t]].mean(1) - c0).reindex(Sg).abs() / fcg.reindex(Sg).abs()).clip(0, None) for t in ['R1W', 'R3W']}
ax = fig.add_subplot(gs[0, 1])
xt = np.arange(2); w = 0.36
for k, (mn, gs) in enumerate(mm.items()):
    ens = [e for e in Sg if es.get(e, '') in set(gs)]
    vals = [qg['R1W'].reindex(ens).mean(), qg['R3W'].reindex(ens).mean()]
    ax.bar(xt + (k - 0.5) * w, vals, w, color=(CI if mn == 'Macrophage' else CE), label=mn)
ax.set_xticks(xt); ax.set_xticklabels(["1 wk", "3 wk"], fontsize=5.6)
ax.set_xlabel("weeks after withdrawal"); ax.set_ylabel("residual $q$ (GSE135462)")
ax.set_ylim(0, 1.0)
ax.legend(loc='upper right', fontsize=5.4, handlelength=1.2, labelspacing=0.3)
S.plab(ax, "b")
S.save(fig, IMG, "FigS3")
print("supp done")
