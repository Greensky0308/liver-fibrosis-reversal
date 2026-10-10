#!/usr/bin/env python3
"""11_compartment_within_mode.py — 绕开"弱药"混淆: 检验"免疫残余 > ECM 残余"
这一分区不对称是否在【每种逆转模式内部】都成立(而非跨模式比谁撤得多).

三种模式 (纤维化=模型vs对照; 逆转=干预vs对照):
  A 毒素撤除  GSE199392  P6 vs P0 ; R6(撤除6w)
  B 饮食撤除  GSE196908 Exp2  NASH-veh vs Chow ; Chow-Reversal
  C 药物治疗  GSE196908 Exp1  NASH-veh12w vs Chow+veh ; Lanifibranor12w / Semaglutide12w
基因集: MSigDB 小鼠 无偏库(免疫/炎症 vs ECM/胶原).
指标: q=|FC_逆转|/|FC_纤维化| (0=完全回落); 比较 各集 mean q.
若 免疫q > ECMq 在 A/B/C 内均成立 => 分区不对称与"药强药弱"无关.
"""
import os, re, gzip
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT  = os.path.join(ROOT, "data")
FC = 1.0

def load_gmt(fn):
    s = {}
    for line in open(os.path.join(DAT, fn)):
        p = line.rstrip("\n").split("\t"); s[p[0]] = [g for g in p[2:] if g]
    return s
mh = load_gmt("mh.all.v2023.2.Mm.symbols.gmt"); m5 = load_gmt("m5.go.bp.v2023.2.Mm.symbols.gmt")
IMMUNE = {
    'INFLAMMATORY_RESPONSE': mh['HALLMARK_INFLAMMATORY_RESPONSE'],
    'TNFA_NFKB': mh['HALLMARK_TNFA_SIGNALING_VIA_NFKB'],
    'IL6_JAK_STAT3': mh['HALLMARK_IL6_JAK_STAT3_SIGNALING'],
    'COMPLEMENT': mh['HALLMARK_COMPLEMENT'],
    'INNATE_IMMUNE': m5['GOBP_INNATE_IMMUNE_RESPONSE'],
    'INFLAMMATORY_BP': m5['GOBP_INFLAMMATORY_RESPONSE'],
    'MACROPHAGE_ACTIVATION': m5['GOBP_MACROPHAGE_ACTIVATION'],
    'LEUKOCYTE_MIGRATION': m5['GOBP_LEUKOCYTE_MIGRATION'],
}
ECM = {
    'EMT': mh['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'],
    'COLLAGEN_FIBRIL_ORG': m5['GOBP_COLLAGEN_FIBRIL_ORGANIZATION'],
    'ECM_ASSEMBLY': m5['GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY'],
    'COLLAGEN_BIOSYNTH': m5['GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS'],
}
gi = pd.read_csv(os.path.join(DAT, "Mus_musculus.gene_info.gz"), sep="\t", low_memory=False)
gi.columns = [c.lstrip('#') for c in gi.columns]
sym2ens = {}
for s, xr in zip(gi['Symbol'], gi['dbXrefs'].fillna('')):
    for m in re.findall(r'Ensembl:(ENSMUSG\d+)', xr): sym2ens.setdefault(s, m)

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
id2g = sm_id2group(os.path.join(DAT, "GSE196908_series_matrix.txt.gz"))
x1 = log2(pd.read_excel(os.path.join(DAT, "GSE199392_processed_gene_expression_matrix.xlsx")).set_index("#ID").astype(float))
def load_exp(f):
    d = pd.read_csv(os.path.join(DAT, f), sep="\t", index_col=0).astype(float)
    d.columns = [id2g.get(str(c), 'NA') for c in d.columns]; return log2(d)
e1 = load_exp("GSE196908_rpkm.txt.gz"); e2 = load_exp("GSE196908_Exp2_rpkm.txt.gz")

def qvec(df, fibcols, ctrlcols, rescols):
    c0 = df[ctrlcols].mean(1); fc = df[fibcols].mean(1) - c0
    S = fc[fc.abs() >= FC].index
    q = ((df[rescols].mean(1) - c0).reindex(S).abs() / fc.reindex(S).abs()).clip(0, None)
    return q

MODES = {
 'A 毒素撤除 (GSE199392 CCl4)': qvec(x1, ['P6-%d'%i for i in [1,2,3,5,6]], ['P0-%d'%i for i in range(1,6)], ['R6-%d'%i for i in range(1,6)]),
 'B 饮食撤除 (Exp2)':           qvec(e2, ['DIO-NASH vehicle w12'], ['Chow vehicle w12'], ['DIO-NASH Chow Reversal w12']),
 'C 药物 Lanifibranor (Exp1)':  qvec(e1, ['DIO-NASH Vehicle 12w'], ['Chow+Vehicle 12w'], ['Lanifibranor 12w']),
 'C 药物 Semaglutide (Exp1)':   qvec(e1, ['DIO-NASH Vehicle 12w'], ['Chow+Vehicle 12w'], ['Semaglutide 12w']),
}

def setq(q, genes):
    hit = [sym2ens[g] for g in genes if g in sym2ens and sym2ens[g] in q.index]
    return q.reindex(hit).dropna()

print(f"{'模式':30s} {'签名n':>6s} | {'免疫集 mean q':>13s} | {'ECM集 mean q':>13s} | {'免疫-ECM':>9s} | {'p(MWU)':>9s}")
for name, q in MODES.items():
    imm = pd.concat([setq(q, g) for g in IMMUNE.values()])
    ecm = pd.concat([setq(q, g) for g in ECM.values()])
    try: _, p = mannwhitneyu(imm, ecm, alternative='greater')
    except ValueError: p = np.nan
    print(f"{name:30s} {len(q):6d} | {imm.mean():13.3f} | {ecm.mean():13.3f} | "
          f"{imm.mean()-ecm.mean():+9.3f} | {p:9.2e}")

print("\n--- 逐集明细 (mean q, 每模式) ---")
for name, q in MODES.items():
    print(f"\n[{name}]  整体签名 q 中位={q.median():.3f}")
    for tag, panel in [('IMMUNE', IMMUNE), ('ECM', ECM)]:
        for gs, genes in panel.items():
            v = setq(q, genes)
            if len(v) >= 5:
                print(f"   {tag:6s} {gs:24s} n={len(v):4d}  q={v.mean():.3f}")
print("\nDONE")
