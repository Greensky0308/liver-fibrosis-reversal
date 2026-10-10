#!/usr/bin/env python3
"""01_reversal_compare.py — 三种肝纤维化逆转模式并排比较 (v2, 修正样本映射)
映射依据 series_matrix 的 !Sample_description(样本数字ID) -> !Sample_characteristics(treatment)
"""
import os, gzip
import numpy as np, pandas as pd
from scipy.stats import ttest_ind
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT=os.path.join(ROOT,"data"); CACHE=os.path.join(ROOT,"out")
os.makedirs(CACHE,exist_ok=True)
FC=1.0; P=0.05
def log2(x): return np.log2(x+1)

def sm_id2group(path):
    import re
    desc=tit=None
    for line in gzip.open(path,'rt',errors='ignore'):
        f=line.rstrip('\n').split('\t')
        if line.startswith('!Sample_description'): desc=f[1:]
        if line.startswith('!Sample_title'): tit=f[1:]
    desc=[d.strip('"') for d in desc]
    def norm(t):
        t=t.strip('"')
        t=re.sub(r'\s*Rep\s*\d+\s*$','',t); t=re.sub(r'^Exp\s*2\s*','',t)
        return t.strip()
    tit=[norm(t) for t in tit]
    return dict(zip(desc,tit))

# ---------- GSE199392 ----------
x1=log2(pd.read_excel(os.path.join(DAT,"GSE199392_processed_gene_expression_matrix.xlsx")).set_index("#ID").astype(float))
g1={'control':['P0-%d'%i for i in range(1,6)],'prog6W':['P6-%d'%i for i in [1,2,3,5,6]],
    'R1':['R1-%d'%i for i in range(1,6)],'R3':['R3-%d'%i for i in range(1,6)],'R6':['R6-%d'%i for i in range(1,6)]}
print("GSE199392:",x1.shape)

# ---------- GSE196908 (both experiments) ----------
id2g=sm_id2group(os.path.join(DAT,"GSE196908_series_matrix.txt.gz"))
def load_exp(f):
    d=pd.read_csv(os.path.join(DAT,f),sep="\t",index_col=0).astype(float)
    g=[id2g.get(str(c),'NA') for c in d.columns]
    d.columns=g
    return log2(d), pd.Series(g)
e1,gE1=load_exp("GSE196908_rpkm.txt.gz")
e2,gE2=load_exp("GSE196908_Exp2_rpkm.txt.gz")
print("Exp1:",e1.shape,"groups:",gE1.value_counts().to_dict())
print("Exp2:",e2.shape,"groups:",gE2.value_counts().to_dict())

def de(df,a,b):
    A=df[a].values; B=df[b].values
    fc=A.mean(1)-B.mean(1); pv=ttest_ind(A,B,axis=1,equal_var=False).pvalue
    return pd.Series(fc,index=df.index), pd.Series(pv,index=df.index)

# ===== 模式A: 毒素撤除 =====
print("\n=== 模式A: 毒素撤除 (GSE199392) ===")
fc,pv=de(x1,g1['prog6W'],g1['control']); sigA=fc[fc.abs()>=FC].index
print(f"纤维化签名(prog6W vs control): {len(sigA)}")
for tp in ['R1','R3','R6']:
    rfc,_=de(x1,g1[tp],g1['control'])
    rev=sigA[rfc.reindex(sigA).abs()<FC/2]
    print(f"  {tp} vs control: 回落 {len(rev)}/{len(sigA)} ({len(rev)/len(sigA):.1%})")

# ===== 模式B: 饮食撤除 =====
print("\n=== 模式B: 饮食撤除 (Exp2) ===")
fc,pv=de(e2,['DIO-NASH vehicle w12'],['Chow vehicle w12']); sigB=fc[fc.abs()>=FC].index
print(f"纤维化签名(NASHveh vs chow): {len(sigB)}")
rfc,_=de(e2,['DIO-NASH Chow Reversal w12'],['Chow vehicle w12'])
revB=sigB[rfc.reindex(sigB).abs()<FC/2]
print(f"  饮食撤除后回落: {len(revB)}/{len(sigB)} ({len(revB)/len(sigB):.1%})")

# ===== 模式C: 药物治疗 =====
print("\n=== 模式C: 药物治疗 (Exp1) ===")
fc,pv=de(e1,['DIO-NASH Vehicle 12w'],['Chow+Vehicle 12w']); sigC=fc[fc.abs()>=FC].index
print(f"纤维化签名(NASHveh12 vs chow): {len(sigC)}")
for tr in ['Lanifibranor 12w','Semaglutide 12w','Lanifibranor 8w','Semaglutide 8w']:
    rfc,_=de(e1,[tr],['Chow+Vehicle 12w'])
    rev=sigC[rfc.reindex(sigC).abs()<FC/2]
    print(f"  {tr} vs chow: 回落 {len(rev)}/{len(sigC)} ({len(rev)/len(sigC):.1%})")

# ===== 跨模式重叠 =====
print("\n=== 纤维化签名重叠 ===")
for a,b,na,nb in [(sigB,sigC,'B','C'),(sigA,sigB,'A','B'),(sigA,sigC,'A','C')]:
    i=len(set(a)&set(b)); print(f"  {na}∩{nb} = {i} (J={i/max(1,len(set(a)|set(b))):.3f})")

# ===== sanity: 修正映射后各组是否自聚 =====
print("\n=== sanity: Exp1 各处理组 within vs between 相关 ===")
v=e1.var(1).sort_values(ascending=False).head(2000).index
C=np.corrcoef(e1.loc[v].values.T); lab=np.array(gE1)
for g in gE1.value_counts().index:
    idx=np.where(lab==g)[0]; oth=np.where(lab!=g)[0]
    if len(idx)<2: continue
    w=np.mean([C[i,j] for i in idx for j in idx if i<j]); b=np.mean([C[i,j] for i in idx for j in oth])
    print(f"  {g:26s} n={len(idx):2d} within={w:.3f} between={b:.3f} gap={w-b:+.3f}")
print("\nDONE")
