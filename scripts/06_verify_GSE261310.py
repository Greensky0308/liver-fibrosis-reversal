#!/usr/bin/env python3
"""06_verify_GSE261310.py — 独立验证: 逆转是否非单调(先恢复再恶化)?
GSE261310: TAA 8w 造模 -> 撤除 1/2/4 周。看全转录组失调度随撤除时间的走势。
对比 GSE199392 的"R3 最优、R6 反弹"。
"""
import os,re
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT=os.path.join(ROOT,"data"); CACHE=os.path.join(ROOT,"out")
d=pd.read_excel(os.path.join(DAT,"GSE261310_Matrix.GRCm38.xlsx"))
d=d.rename(columns={d.columns[0]:"sym",d.columns[1]:"ens"})
keep=["sym","ens"]+[c for c in d.columns if re.match(r'^(WT|T8|T12|T8F)',str(c))]
d=d[keep].set_index("ens")
symser=d["sym"]
d=d.drop(columns=["sym"])
X=np.log2(d.astype(float)+1)
def cols(p): return [c for c in X.columns if re.match(p,c)]
G={'WT':cols(r'WT_'),'T8':cols(r'T8_\d'),'T12':cols(r'T12_'),
   'F1':cols(r'T8F1_\d'),'F2':cols(r'T8F2_\d'),'F4':cols(r'T8F4_\d'),
   'F1C':cols(r'T8F1C'),'F2C':cols(r'T8F2C'),'F4C':cols(r'T8F4C')}
print({k:len(v) for k,v in G.items()})
c0=X[G['WT']].mean(1)
print("\n=== 全转录组失调度 (vs WT) ===")
for k in ['T8','T12','F1','F2','F4']:
    dmean=X[G[k]].mean(1)-c0
    print(f"  {k:4s} mean|log2FC|={dmean.abs().mean():.3f} | |FC|>1 基因数={(dmean.abs()>1).sum()}")

# 纤维化签名 (T8 vs WT), 看撤除后走势
sig=(X[G['T8']].mean(1)-c0).abs(); sig=sig[sig>1].index
print(f"\nT8 纤维化签名: {len(sig)}")
res={}
for k in ['T8','T12','F1','F2','F4']:
    dv=(X[G[k]].mean(1)-c0).reindex(sig)
    q=(dv.abs()/(X[G['T8']].mean(1)-c0).reindex(sig).abs()).clip(0,1)
    res[k]=q.median()
print("签名中位分辨率 q (1=仍纤维化,0=完全回落):")
for k,v in res.items(): print(f"  {k:4s} q_median={v:.2f}")

# 关键标志物
mk=['Col1a1','Col3a1','Acta2','Timp1','Postn','Lox','Cd36','Cd163','Mmp7']
s2e={v:k for k,v in symser.items()}
print("\n关键标志物 (log2):")
for g in mk:
    e=s2e.get(g)
    if e in X.index:
        vals=[X.loc[e,G[k]].mean() for k in ['WT','T8','F1','F2','F4']]
        print(f"  {g:8s} WT={vals[0]:.1f} T8={vals[1]:.1f} F1={vals[2]:.1f} F2={vals[3]:.1f} F4={vals[4]:.1f}")
print("\nDONE")
