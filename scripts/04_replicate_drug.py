#!/usr/bin/env python3
"""04_replicate_drug.py — 第二个药物逆转臂 (GSE272274, mARC1-siRNA)
检验: 共识签名里 "cause_only" 基因(撤因可撤、lanifibranor 撤不掉) 是否被另一个药也撤不掉?
若两个独立药物都撤不掉 -> 不对称更可能真实, 而非单一药太弱。
"""
import os, re, gzip
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DAT=os.path.join(ROOT,"data"); CACHE=os.path.join(ROOT,"out")
FC=1.0
def log2(x): return np.log2(x+1)

d=pd.read_csv(os.path.join(DAT,"GSE272274_RPKMCodingGenes.xls.gz"),sep="\t")
gc=d.columns[-1]; d=d.rename(columns={gc:"ens"}).set_index("ens")
def grp(c):
    c=str(c).replace("C57BL6J_","").replace("_FPKM","")
    return re.sub(r"_S\d+$","",c)
d.columns=[grp(c) for c in d.columns]
print("GSE272274 groups:", pd.Series(d.columns).value_counts().to_dict())
X=log2(d.astype(float))
gchow=[c for c in X.columns if "Chow_Diet" in c]
gveh =[c for c in X.columns if "PBS" in c and "Chow" not in c]
gdrg =[c for c in X.columns if "siMARC1" in c]
print(f"chow={len(gchow)} veh={len(gveh)} drug={len(gdrg)}")
def de(a,b): return pd.Series(X[a].values.mean(1)-X[b].values.mean(1),index=X.index)
fc_fib=de(gveh,gchow); fc_res=de(gdrg,gchow)
sig=fc_fib[fc_fib.abs()>=FC].index
print(f"GSE272274 纤维化签名(HFD vs chow): {len(sig)}")
r=(1-(fc_res.reindex(sig).abs()/fc_fib.reindex(sig).abs())).clip(0,1)
print(f"mARC1-siRNA 可逆比中位: {r.median():.2f} | reversed(>=.7): {(r>=.7).mean():.1%}")

# 符号
gi=pd.read_csv(os.path.join(DAT,"Mus_musculus.gene_info.gz"),sep="\t",low_memory=False,header=0)
gi.columns=[c.lstrip('#') for c in gi.columns]
sym={}
for sm,xr in zip(gi['Symbol'],gi['dbXrefs'].fillna('')):
    for m in re.findall(r'Ensembl:(ENSMUSG\d+)',xr): sym[m]=sm
r.index=[sym.get(g,g) for g in r.index]

# 取共识 cause_only / persistent_all 基因, 看新药能否撤
D=pd.read_csv(os.path.join(CACHE,"consensus_reversibility.csv"),index_col=0)
for cls in ['cause_only','persistent_all','reversed_by_both']:
    genes=list(D[D['class']==cls]['symbol'])
    genes=[g for g in genes if g in r.index]
    if not genes: print(f"{cls}: 无可比基因"); continue
    rr=r.reindex(genes)
    print(f"\n{cls} (n={len(genes)}) 在 mARC1-siRNA 下的可逆比: 中位 {rr.median():.2f} | reversed(>=.7) {(rr>=.7).mean():.1%}")
    print("  ", dict(zip(genes, np.round(rr.values,2))))
print("\nDONE")
