#!/usr/bin/env python3
"""14_composition_vs_state.py — 判断"免疫残余"= 细胞变多(构成) 还是 细胞状态没关.

① marker-based 细胞丰度轨迹 (P0/P6/R1/R3/R6): 免疫细胞在撤除后回落了吗?
② 免疫模块拆: 身份标记基因(随丰度) vs 效应基因(随状态)
③ 用丰度分数回归校正每个签名基因 -> 重算残余 q -> 免疫 vs ECM 对比是否仍成立
数据: GSE199392 (CCl4, 撤除 1/3/6 周). 复用 GSE135462 交叉验证.
"""
import os, re
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT  = os.path.join(ROOT, "data")
FC = 1.0

gi = pd.read_csv(os.path.join(DAT, "Mus_musculus.gene_info.gz"), sep="\t", low_memory=False)
gi.columns = [c.lstrip('#') for c in gi.columns]
sym2ens = {}
for s, xr in zip(gi['Symbol'], gi['dbXrefs'].fillna('')):
    for m in re.findall(r'Ensembl:(ENSMUSG\d+)', xr): sym2ens.setdefault(s, m)

# ---- 细胞类型 marker (肝) ----
CELL = {
 'Mac_Kupffer': ['Cd68','Csf1r','Adgre1','Mrc1','Cd163','Lyz2','Msr1','Trem2','Vsig4','Marco','Timd4','Clec4f','Fcgr1'],
 'Monocyte':    ['Ccr2','Ly6c2','Chil3','S100a8','S100a9','Fcnb','Ms4a6c'],
 'T_cell':      ['Cd3e','Cd3d','Cd3g','Cd8a','Cd4','Trbc2','Lck','Thy1'],
 'NK':          ['Nkg7','Klrb1c','Ncr1','Klrd1','Prf1'],
 'B_cell':      ['Cd79a','Cd79b','Ms4a1','Cd19','Ighm'],
 'Neutrophil':  ['Ly6g','Mpo','Elane','Retnlg','S100a8','S100a9'],
 'HSC_mesench': ['Acta2','Des','Pdgfrb','Gfap','Gli1','Col1a1','Col3a1'],
 'Hepatocyte':  ['Alb','Apoa1','Apoa2','Ttr','Serpina1a','Cyp2e1','Asgr1','Fga'],
 'Endothelial': ['Pecam1','Cdh5','Kdr','Stab2','Clec4g','Lyve1'],
}
IMMUNE_CELLS = ['Mac_Kupffer','Monocyte','T_cell','NK','B_cell','Neutrophil']

# ---- 无偏基因集: 免疫/炎症 vs ECM ----
def load_gmt(fn):
    s = {}
    for line in open(os.path.join(DAT, fn)):
        p = line.rstrip("\n").split("\t"); s[p[0]] = [g for g in p[2:] if g]
    return s
mh = load_gmt("mh.all.v2023.2.Mm.symbols.gmt"); m5 = load_gmt("m5.go.bp.v2023.2.Mm.symbols.gmt")
IMMUNE_SET = {'INFLAMM': mh['HALLMARK_INFLAMMATORY_RESPONSE'], 'TNFA': mh['HALLMARK_TNFA_SIGNALING_VIA_NFKB'],
              'IL6': mh['HALLMARK_IL6_JAK_STAT3_SIGNALING'], 'COMPL': mh['HALLMARK_COMPLEMENT'],
              'INNATE': m5['GOBP_INNATE_IMMUNE_RESPONSE'], 'MAC_ACT': m5['GOBP_MACROPHAGE_ACTIVATION'],
              'LEUK_MIG': m5['GOBP_LEUKOCYTE_MIGRATION']}
ECM_SET = {'EMT': mh['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'],
           'COL_FIBRIL': m5['GOBP_COLLAGEN_FIBRIL_ORGANIZATION'],
           'ECM_ASM': m5['GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY'],
           'COL_BIOSYN': m5['GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS']}
# 免疫"身份标记"(定义丰度的基因) 全集, 用于把模块拆成 身份 vs 效应
IDENTITY = set(sum((CELL[c] for c in IMMUNE_CELLS), []))

def analyze(name, X, G, times, last, ctrl, fib):
    c0 = X[G[ctrl]].mean(1)
    # ---- ① 细胞丰度轨迹 ----
    print(f"\n{'='*74}\n### {name}  细胞类型分数轨迹 (marker 均值, 高=该细胞多)")
    Xsym = {s: X.loc[e] for s, e in sym2ens.items() if e in X.index}
    scores = {}
    for ct, gs in CELL.items():
        vs = [Xsym[g] for g in gs if g in Xsym]
        if vs: scores[ct] = pd.concat(vs, axis=1).mean(1)
    sc = pd.DataFrame(scores)              # 样本 x 细胞类型
    tr = pd.DataFrame({t: sc.loc[G[t]].mean() for t in [ctrl, fib] + times})
    print(tr.round(2).to_string())
    imm_sc = sc[IMMUNE_CELLS].sum(1)

    # ---- ② 签名 + 免疫模块拆 身份 vs 效应 ----
    fc = X[G[fib]].mean(1) - c0
    S = fc[fc.abs() >= FC].index
    ens2sym = {e: next((s for s, e2 in sym2ens.items() if e2 == e), '') for e in S}

    def qvals(cols):
        return ((X[G[cols]].mean(1) - c0).reindex(S).abs() / fc.reindex(S).abs()).clip(0, None)

    def setg(sets):
        out = []
        for gs in sets.values():
            out += [e for e in S if ens2sym.get(e, '') in set(gs)]
        return sorted(set(out))

    imm_g = setg(IMMUNE_SET); ecm_g = setg(ECM_SET)
    imm_ident = [e for e in imm_g if ens2sym.get(e, '') in IDENTITY]
    imm_eff = [e for e in imm_g if ens2sym.get(e, '') not in IDENTITY]
    q = qvals(last)
    print(f"\n  签名 {len(S)} | 免疫集 {len(imm_g)}(身份{len(imm_ident)}/效应{len(imm_eff)}) | ECM集 {len(ecm_g)}")
    print(f"  末点({last}) 残余 q:  免疫全部={q.reindex(imm_g).mean():.3f}  "
          f"身份={q.reindex(imm_ident).mean():.3f}  效应={q.reindex(imm_eff).mean():.3f}  ECM={q.reindex(ecm_g).mean():.3f}")
    print(f"  MWU 免疫全部>ECM p={mannwhitneyu(q.reindex(imm_g), q.reindex(ecm_g), alternative='greater').pvalue:.2e}")

    # ---- ③ 丰度校正 ----
    # 每个签名基因 ~ 免疫丰度分数 (全体样本), 取残差 -> 重算 q
    allc = G[ctrl] + G[fib] + sum((G[t] for t in times), [])
    A = imm_sc.reindex(allc).values
    A = (A - A.mean()) / (A.std() + 1e-9)
    A = np.column_stack([np.ones(len(A)), A])
    Xs = X.reindex(S)[allc].values        # genes x samples
    beta, *_ = np.linalg.lstsq(A, Xs.T, rcond=None)
    resid = (Xs.T - A @ beta).T           # 校正后
    Xr = pd.DataFrame(resid, index=S, columns=allc)
    c0r = Xr[G[ctrl]].mean(1); fcr = Xr[G[fib]].mean(1) - c0r
    qr = ((Xr[G[last]].mean(1) - c0r).abs() / fcr.abs()).clip(0, None)
    print(f"\n  [丰度校正后] 免疫全部={qr.reindex(imm_g).mean():.3f} 身份={qr.reindex(imm_ident).mean():.3f} "
          f"效应={qr.reindex(imm_eff).mean():.3f} ECM={qr.reindex(ecm_g).mean():.3f}")
    p = mannwhitneyu(qr.reindex(imm_g), qr.reindex(ecm_g), alternative='greater').pvalue
    print(f"  [丰度校正后] MWU 免疫>ECM p={p:.2e}")
    return q, qr, imm_g, ecm_g, imm_ident, imm_eff

# GSE199392
X1 = np.log2(pd.read_excel(os.path.join(DAT, "GSE199392_processed_gene_expression_matrix.xlsx")).set_index("#ID").astype(float) + 1)
G1 = {'ctrl': ['P0-%d'%i for i in range(1,6)], 'fib': ['P6-%d'%i for i in [1,2,3,5,6]],
      'R1': ['R1-%d'%i for i in range(1,6)], 'R3': ['R3-%d'%i for i in range(1,6)], 'R6': ['R6-%d'%i for i in range(1,6)]}
analyze("GSE199392 CCl4 (撤除 1/3/6 周)", X1, G1, ['R1','R3','R6'], 'R6', 'ctrl', 'fib')

# GSE135462
dg = pd.read_csv(os.path.join(DAT, "GSE135462_FPKM_recovered.csv")).set_index("ens")
X2 = np.log2(dg[[f"S{i}" for i in range(1,14)]].astype(float) + 1)
G2 = {'ctrl': ['S11','S12','S13'], 'fib': ['S5','S6'], 'R1W': ['S7','S8'], 'R3W': ['S9','S10']}
analyze("GSE135462 CCl4 (撤除 1/3 周)", X2, G2, ['R1W','R3W'], 'R3W', 'ctrl', 'fib')
print("\nDONE")
