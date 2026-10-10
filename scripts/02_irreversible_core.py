#!/usr/bin/env python3
"""02_irreversible_core.py — 不可逆核心 + 模式特异逆转
对每个模型：纤维化诱导 FC_fib，逆转后残余 FC_res；可逆比 r = 1-|FC_res|/|FC_fib|
分类: reversed(r>=.7) / partial(.3-.7) / persistent(<.3)
跨模式: 不可逆核心 = 在 A/B/C 都被诱导且都 persistent 的基因
"""
import os, gzip, re
import numpy as np, pandas as pd
from scipy.stats import ttest_ind

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT=os.path.join(ROOT,"data"); CACHE=os.path.join(ROOT,"out")
FC=1.0; P=0.05
def log2(x): return np.log2(x+1)
def sm_id2group(path):
    desc=tit=None
    for line in gzip.open(path,'rt',errors='ignore'):
        f=line.rstrip('\n').split('\t')
        if line.startswith('!Sample_description'): desc=f[1:]
        if line.startswith('!Sample_title'): tit=f[1:]
    desc=[d.strip('"') for d in desc]
    tit=[re.sub(r'^Exp\s*2\s*','',re.sub(r'\s*Rep\s*\d+\s*$','',t.strip('"'))) for t in tit]
    return dict(zip(desc,tit))
id2g=sm_id2group(os.path.join(DAT,"GSE196908_series_matrix.txt.gz"))
def load_exp(f):
    d=pd.read_csv(os.path.join(DAT,f),sep="\t",index_col=0).astype(float)
    d.columns=[id2g.get(str(c),'NA') for c in d.columns]; return log2(d)

x1=log2(pd.read_excel(os.path.join(DAT,"GSE199392_processed_gene_expression_matrix.xlsx")).set_index("#ID").astype(float))
e1=load_exp("GSE196908_rpkm.txt.gz"); e2=load_exp("GSE196908_Exp2_rpkm.txt.gz")

def de(df,a,b):
    A=df[a].values; B=df[b].values
    return pd.Series(A.mean(1)-B.mean(1),index=df.index)

MODELS={
 'A_toxinW':  (x1, ['P6-%d'%i for i in [1,2,3,5,6]], ['P0-%d'%i for i in range(1,6)], ['R6-%d'%i for i in range(1,6)]),
 'B_dietW':   (e2, ['DIO-NASH vehicle w12'], ['Chow vehicle w12'], ['DIO-NASH Chow Reversal w12']),
 'C_drug':    (e1, ['DIO-NASH Vehicle 12w'], ['Chow+Vehicle 12w'], ['Lanifibranor 12w']),
}
res={}
for name,(df,fib,ctl,rev) in MODELS.items():
    fc_fib=de(df,fib,ctl); fc_res=de(df,rev,ctl)
    sig=fc_fib[fc_fib.abs()>=FC].index
    r=(1 - (fc_res.reindex(sig).abs()/fc_fib.reindex(sig).abs())).clip(0,1)
    res[name]=r
    cls=pd.cut(r,[-0.01,0.3,0.7,1.01],labels=['persistent','partial','reversed'])
    print(f"{name}: signature={len(sig)} | 可逆比中位={r.median():.2f} | {cls.value_counts().to_dict()}")

# 跨模式 (三数据集共基因)
common=set(x1.index)&set(e1.index)&set(e2.index)
print("\n三数据集共基因:",len(common))
ra=res['A_toxinW'].reindex(sorted(common)); rb=res['B_dietW'].reindex(sorted(common)); rc=res['C_drug'].reindex(sorted(common))
# 不可逆核心: 三模式都存在(非NaN) 且 都 persistent(<0.3)
core=ra.notna()&rb.notna()&rc.notna()&(ra<0.3)&(rb<0.3)&(rc<0.3)
# 逆转核心: 至少两模式 reversed(>=0.7)
revc=((ra>=0.7).astype(int)+(rb>=0.7).astype(int)+(rc>=0.7).astype(int))>=2
print(f"\n不可逆核心 (三模式都 persistent<0.3): {core.sum()}")
print(f"可逆核心 (>=2 模式 reversed>=0.7): {revc.sum()}")
# 模式特异: 撤除可逆但药物不可
withdr_rev_drug_pers = ((ra>=0.7)|(rb>=0.7)) & (rc<0.3) & ra.notna()&rb.notna()&rc.notna()
drug_rev_withdr_pers = (rc>=0.7) & ((ra<0.3)&(rb<0.3)) & ra.notna()&rb.notna()&rc.notna()
print(f"撤除可逆 & 药物撤不掉: {withdr_rev_drug_pers.sum()}")
print(f"药物可逆 & 撤除撤不掉: {drug_rev_withdr_pers.sum()}")

# 符号映射 (gene_info dbXrefs)
sym={}
gi=pd.read_csv(os.path.join(DAT,"Mus_musculus.gene_info.gz"),sep="\t",low_memory=False,header=0)
gi.columns=[c.lstrip('#') for c in gi.columns]
for sm,xr in zip(gi['Symbol'],gi['dbXrefs'].fillna('')):
    for m in re.findall(r'Ensembl:(ENSMUSG\d+)',xr): sym[m]=sm
print("symbol mapped:",len(sym))
def show(idx):
    return [sym.get(g,g) for g in idx]
print("\n不可逆核心 基因(前30):", show(ra.index[core][:30]))
print("\n撤除可逆&药物撤不掉 (前30):", show(ra.index[withdr_rev_drug_pers][:30]))
print("\n药物可逆&撤除撤不掉 (前30):", show(ra.index[drug_rev_withdr_pers][:30]))
pd.DataFrame({'A_toxinW_r':ra,'B_dietW_r':rb,'C_drug_r':rc}).to_csv(os.path.join(CACHE,"reversibility_ratio.csv"))
print("\nDONE")
