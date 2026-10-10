#!/usr/bin/env python3
"""03_consensus_clusters.py — 共识纤维化签名 + 可逆性聚类 + 通路注释
共识签名 = 在 >=2 个模型中 |FC_fib|>=1 且同向的基因(三数据集共有)。
对共识基因算三模式可逆比 r，聚类为:
  reversed_by_all / cause_only(撤除可逆&药物不可) / drug_only / persistent_all / partial
并按 ECM/巨噬/HSC 模块注释。
"""
import os, gzip, re
import numpy as np, pandas as pd
from scipy.stats import ttest_ind

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT=os.path.join(ROOT,"data"); CACHE=os.path.join(ROOT,"out")
FC=1.0
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
    A=df[a].values; B=df[b].values; return pd.Series(A.mean(1)-B.mean(1),index=df.index)
M={
 'A_toxinW':(x1,['P6-%d'%i for i in [1,2,3,5,6]],['P0-%d'%i for i in range(1,6)],['R6-%d'%i for i in range(1,6)]),
 'B_dietW': (e2,['DIO-NASH vehicle w12'],['Chow vehicle w12'],['DIO-NASH Chow Reversal w12']),
 'C_drug':  (e1,['DIO-NASH Vehicle 12w'],['Chow+Vehicle 12w'],['Lanifibranor 12w']),
}
fib={};resd={}
for k,(df,f,c,r) in M.items():
    fib[k]=de(df,f,c); resd[k]=de(df,r,c)

# 共识签名: >=2 模型 |FC|>=1 且同向
common=set(x1.index)&set(e1.index)&set(e2.index)
sig=[]
for g in common:
    vals=[fib[k].get(g,np.nan) for k in M]
    up=sum(1 for v in vals if v==v and v>=FC); dn=sum(1 for v in vals if v==v and v<=-FC)
    if up>=2: sig.append(g)
    elif dn>=2: sig.append(g)
print(f"共识纤维化签名 (>=2模型,|FC|>=1,同向): {len(sig)}")

def ratio(k):
    f=fib[k].reindex(sig); r=resd[k].reindex(sig)
    return (1-(r.abs()/f.abs())).clip(0,1)
rA,rB,rC=ratio('A_toxinW'),ratio('B_dietW'),ratio('C_drug')
D=pd.DataFrame({'rA_toxinW':rA,'rB_dietW':rB,'rC_drug':rC,
                'fibA':fib['A_toxinW'].reindex(sig),'fibB':fib['B_dietW'].reindex(sig),'fibC':fib['C_drug'].reindex(sig)})
cause=np.maximum(rA,rB)
D['cause_r']=cause; D['drug_r']=rC
D.dropna(subset=['rA_toxinW','rB_dietW','rC_drug'],how='all',inplace=True)

def classify(row):
    ca,dr=row['cause_r'],row['drug_r']
    if pd.isna(ca) or pd.isna(dr): return 'na'
    if ca>=0.7 and dr>=0.7: return 'reversed_by_both'
    if ca>=0.7 and dr<0.3:  return 'cause_only'
    if dr>=0.7 and ca<0.3:  return 'drug_only'
    if ca<0.3 and dr<0.3:   return 'persistent_all'
    return 'partial'
D['class']=D.apply(classify,axis=1)
print("\n聚类分布:", D['class'].value_counts().to_dict())

# 符号
gi=pd.read_csv(os.path.join(DAT,"Mus_musculus.gene_info.gz"),sep="\t",low_memory=False,header=0)
gi.columns=[c.lstrip('#') for c in gi.columns]
sym={}
for sm,xr in zip(gi['Symbol'],gi['dbXrefs'].fillna('')):
    for m in re.findall(r'Ensembl:(ENSMUSG\d+)',xr): sym[m]=sm
D['symbol']=[sym.get(g,'') for g in D.index]
D.to_csv(os.path.join(CACHE,"consensus_reversibility.csv"))

# 模块注释
MOD={
 'ECM_structural':['Col1a1','Col1a2','Col3a1','Col4a1','Col4a2','Col5a1','Col5a2','Col6a1','Col6a2','Col6a3','Col8a1','Eln','Fn1','Lama4','Lamb1','Lamb2','Lamc1','Vcan','Fmod','Bgn','Dcn','Sparc','Postn'],
 'MMP_TIMP':['Mmp2','Mmp7','Mmp9','Mmp12','Mmp13','Mmp14','Mmp19','Timp1','Timp2','Timp3'],
 'crosslink_LOX':['Lox','Loxl1','Loxl2','Loxl3','Loxl4'],
 'HSC_activation':['Acta2','Des','Pdgfrb','Tgfb1','Tgfbr1','Spp1','Tagln','Col1a1','Timp1','Ctgf','Ccn2'],
 'macrophage':['Cd68','Cd163','Mrc1','Lyve1','Ccl2','Ccl6','Ccl7','Adgre1','Trem2','Cd36','Csf1r','Msr1','Spp1'],
}
allm=set(sum(MOD.values(),[]))
print("\n各模块基因在共识签名中的代表 & 类别:")
for m,gs in MOD.items():
    hit=[g for g in gs if g in set(D['symbol'])]
    print(f"  [{m}] ({len(hit)}/{len(gs)} in sig)")
    for g in hit:
        row=D[D['symbol']==g].iloc[0]
        print(f"     {g:10s} {row['class']:16s} rA={row['rA_toxinW']:.2f} rB={row['rB_dietW']:.2f} rC={row['rC_drug']:.2f}")

print("\n--- persistent_all 前25 ---")
print(D[D['class']=='persistent_all'].sort_values('cause_r')[['symbol','rA_toxinW','rB_dietW','rC_drug']].head(25).to_string(index=False))
print("\n--- cause_only 前25 ---")
print(D[D['class']=='cause_only'].sort_values('drug_r')[['symbol','rA_toxinW','rB_dietW','rC_drug']].head(25).to_string(index=False))
print("\nDONE")
