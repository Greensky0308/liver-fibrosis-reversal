#!/usr/bin/env python3
"""05_kinetics.py — 逆转动力学: 单模型内 (GSE199392) 定义不归点
对每个纤维化基因(P6 vs P0,|FC|>=1): 分辨率 q(t)=|FC_t|/|FC_P6| (1=仍纤维化, 0=完全回落)
分类: fast(R1<=0.5) / late(R3<=0.5) / slow(R6<=0.5) / rebound(触底后 R6 回升>0.3) / persistent(min q>0.5)
"""
import os,re
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT=os.path.join(ROOT,"data"); CACHE=os.path.join(ROOT,"out")
FC=1.0
X=np.log2(pd.read_excel(os.path.join(DAT,"GSE199392_processed_gene_expression_matrix.xlsx")).set_index("#ID").astype(float)+1)
G={'P0':['P0-%d'%i for i in range(1,6)],'P6':['P6-%d'%i for i in [1,2,3,5,6]],
   'R1':['R1-%d'%i for i in range(1,6)],'R3':['R3-%d'%i for i in range(1,6)],'R6':['R6-%d'%i for i in range(1,6)]}
def dev(k): return X[G[k]].mean(1)-X[G['P0']].mean(1)
dP6,dR1,dR3,dR6=dev('P6'),dev('R1'),dev('R3'),dev('R6')
sig=dP6[dP6.abs()>=FC].index
q=pd.DataFrame({'R1':dR1.reindex(sig).abs()/dP6.reindex(sig).abs(),
                'R3':dR3.reindex(sig).abs()/dP6.reindex(sig).abs(),
                'R6':dR6.reindex(sig).abs()/dP6.reindex(sig).abs()}).clip(0,None)
def cls(r):
    if r['R1']<=0.5: base='fast'
    elif r['R3']<=0.5: base='late'
    elif r['R6']<=0.5: base='slow'
    elif min(r['R1'],r['R3'],r['R6'])>0.5: base='persistent'
    else: base='partial'
    if base in ('fast','late','slow') and (r['R6']-min(r['R1'],r['R3'],r['R6']))>0.3: base='rebound'
    return base
q['class']=q.apply(cls,axis=1)
print("纤维化签名:",len(sig))
print("动力学分类:", q['class'].value_counts().to_dict())
# 符号
gi=pd.read_csv(os.path.join(DAT,"Mus_musculus.gene_info.gz"),sep="\t",low_memory=False,header=0); gi.columns=[c.lstrip('#') for c in gi.columns]
sym={}
for s,xr in zip(gi['Symbol'],gi['dbXrefs'].fillna('')):
    for m in re.findall(r'Ensembl:(ENSMUSG\d+)',xr): sym[m]=s
q['symbol']=[sym.get(g,'') for g in q.index]
q.to_csv(os.path.join(CACHE,"kinetics.csv"))
MOD={'ECM':['Col1a1','Col1a2','Col3a1','Col4a1','Col5a1','Col5a2','Col6a1','Col6a2','Postn','Fn1','Bgn','Sparc','Vcan','Fmod','Eln'],
     'MMP_TIMP':['Mmp2','Mmp7','Mmp9','Mmp12','Mmp13','Mmp14','Timp1','Timp2','Timp3'],
     'LOX':['Lox','Loxl1','Loxl2','Loxl3','Loxl4'],
     'HSC':['Acta2','Des','Pdgfrb','Tgfb1','Spp1','Tagln','Ccn2','Ctgf'],
     'MACROPHAGE':['Cd68','Cd163','Mrc1','Lyve1','Ccl2','Ccl6','Ccl7','Adgre1','Trem2','Cd36','Csf1r','Msr1'],
     'LIPID_METAB':['Cd36','Fgf21','Cidea','Cidec','Mogat1','Gpat3','Plin4','G0s2']}
for m,gs in MOD.items():
    hit=[g for g in gs if g in set(q['symbol'])]
    print(f"\n[{m}] ({len(hit)}/{len(gs)})")
    for g in hit:
        row=q[q['symbol']==g].iloc[0]
        print(f"   {g:9s} {row['class']:11s} q: R1={row['R1']:.2f} R3={row['R3']:.2f} R6={row['R6']:.2f}")
print("\n--- persistent (永不回落) top30 by q_min ---")
pers=q[q['class']=='persistent'].copy(); pers['qmin']=pers[['R1','R3','R6']].min(1)
print(pers.sort_values('qmin',ascending=False).head(30)[['symbol','R1','R3','R6']].to_string(index=False))
print("\nDONE")
