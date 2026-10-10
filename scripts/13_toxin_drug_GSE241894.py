#!/usr/bin/env python3
"""13_toxin_drug_GSE241894.py — 解混淆: "毒物模型 + 用药"下残余是否仍是免疫?
GSE241894: CCl4 模型 + 4MU(HAS2 药理抑制) / siHAS2. 健康对照 CNTRL.
若毒物模型"用药"也留免疫残余 => 决定残余的是【损伤类型(毒)】, 而非【撤因】.
对照: 同数据集的撤除式无; 但可与 CCl4 撤除(GSE199392/135462, 免疫残余)对比.
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
IMMUNE = {'INFLAMMATORY_RESPONSE': mh['HALLMARK_INFLAMMATORY_RESPONSE'],
          'TNFA_NFKB': mh['HALLMARK_TNFA_SIGNALING_VIA_NFKB'],
          'IL6_JAK_STAT3': mh['HALLMARK_IL6_JAK_STAT3_SIGNALING'],
          'COMPLEMENT': mh['HALLMARK_COMPLEMENT'],
          'INNATE_IMMUNE': m5['GOBP_INNATE_IMMUNE_RESPONSE'],
          'INFLAMMATORY_BP': m5['GOBP_INFLAMMATORY_RESPONSE'],
          'MACROPHAGE_ACTIVATION': m5['GOBP_MACROPHAGE_ACTIVATION'],
          'LEUKOCYTE_MIGRATION': m5['GOBP_LEUKOCYTE_MIGRATION']}
ECM = {'EMT': mh['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'],
       'COLLAGEN_FIBRIL_ORG': m5['GOBP_COLLAGEN_FIBRIL_ORGANIZATION'],
       'ECM_ASSEMBLY': m5['GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY'],
       'COLLAGEN_BIOSYNTH': m5['GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS']}

gi = pd.read_csv(os.path.join(DAT, "Mus_musculus.gene_info.gz"), sep="\t", low_memory=False)
gi.columns = [c.lstrip('#') for c in gi.columns]
ens2sym = {}
for s, xr in zip(gi['Symbol'], gi['dbXrefs'].fillna('')):
    for m in re.findall(r'Ensembl:(ENSMUSG\d+)', xr): ens2sym[m] = s

d = pd.read_csv(os.path.join(DAT, "GSE241894_counts.csv.gz"), index_col=0)
X = np.log2(d.astype(float) + 1)
def cols(p): return [c for c in X.columns if c.replace(' ','').startswith(p)]
G = {'CNTRL': cols('CNTRL.'), '4MU': cols('4MU.'), 'CCl4': cols('CCl4.'),
     'CCl4_4MU': [c for c in X.columns if 'CCl4+4MU' in c.replace(' ','')],
     'siLuc': cols('siLuc.'), 'CCl4_siLuc': [c for c in X.columns if 'CCl4+siLuc' in c.replace(' ','')],
     'siHAS2': cols('siHAS2.'), 'CCl4_siHAS2': [c for c in X.columns if 'CCl4+siHAS2' in c.replace(' ','')]}
print({k: len(v) for k, v in G.items()})

def setq(q, genes):
    hit = [e for g in genes if (e := next((k for k, v in ens2sym.items() if v == g), None)) and e in q.index]
    return q.reindex(hit).dropna()

# 纤维化签名: 所有 CCl4(模型) vs 所有对照
fibcols = G['CCl4'] + G['CCl4_siLuc']
ctrlcols = G['CNTRL'] + G['siLuc']
c0 = X[ctrlcols].mean(1); fc = X[fibcols].mean(1) - c0
S = fc[fc.abs() >= FC].index
print(f"\n纤维化签名 (CCl4 vs 对照) |FC|>=1: {len(S)}")

ARMS = {'撤除式(对照, 见 GSE199392/135462)': None,
        'CCl4+4MU (药理抑制)': G['CCl4_4MU'],
        'CCl4+siHAS2 (基因抑制)': G['CCl4_siHAS2'],
        'CCl4-only (无治疗)': G['CCl4']}
print(f"\n{'臂':32s} {'q中位':>7s} | {'免疫q':>7s} {'ECMq':>7s} {'免疫-ECM':>9s} {'p(MWU)':>9s}")
for name, cols_ in ARMS.items():
    if name.startswith('撤除'): continue
    q = ((X[cols_].mean(1) - c0).reindex(S).abs() / fc.reindex(S).abs()).clip(0, None)
    imm = pd.concat([setq(q, g) for g in IMMUNE.values()])
    ecm = pd.concat([setq(q, g) for g in ECM.values()])
    p = mannwhitneyu(imm, ecm, alternative='greater').pvalue
    print(f"{name:32s} {q.median():7.3f} | {imm.mean():7.3f} {ecm.mean():7.3f} {imm.mean()-ecm.mean():+9.3f} {p:9.2e}")

print("\n--- 逐集 (CCl4+4MU 臂) ---")
q = ((X[G['CCl4_4MU']].mean(1) - c0).reindex(S).abs() / fc.reindex(S).abs()).clip(0, None)
for tag, panel in [('IMMUNE', IMMUNE), ('ECM', ECM)]:
    for gs, genes in panel.items():
        v = setq(q, genes)
        if len(v) >= 5: print(f"   {tag:6s} {gs:24s} n={len(v):4d} q={v.mean():.3f}")
print("\nDONE")
