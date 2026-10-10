#!/usr/bin/env python3
"""17_singlecell_combat.py — 与 16 同数据, 加 ComBat 批次整合(按样本), 重聚类再算比例.
用于检验"未整合作注释"的细胞比例是否由批次/样本驱动。
"""
import os, re
import numpy as np, pandas as pd
import scipy.io, scipy.sparse as sp
import scanpy as sc, anndata as ad

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SC   = os.path.join(ROOT, "data/GSE295293")
OUT  = os.path.join(ROOT, "out")
sc.settings.verbosity = 1

SAMPLES = {
 'control':    ['GSM8944960_Fixed_RNA_control_wt1198','GSM8944961_Fixed_RNA_control_wt1199','GSM8944962_Fixed_RNA_control_wt1200'],
 'fibrosis':   ['GSM8944963_Fixed_RNA_firosis_wt11','GSM8944964_Fixed_RNA_firosis_wt13','GSM8944965_Fixed_RNA_firosis_wt14'],
 'regression': ['GSM8944966_Fixed_RNA_regression_wt1815','GSM8944967_Fixed_RNA_regression_wt1816','GSM8944968_Fixed_RNA_regression_wt1817'],
}
GROUPS = {s: g for g, ss in SAMPLES.items() for s in ss}

def load(s):
    feat = pd.read_csv(f"{SC}/{s}_features.tsv.gz", sep="\t", header=None)
    bar = pd.read_csv(f"{SC}/{s}_barcodes.tsv.gz", sep="\t", header=None)
    M = scipy.io.mmread(f"{SC}/{s}_matrix.mtx.gz").tocsc().T
    a = ad.AnnData(X=sp.csr_matrix(M), obs=pd.DataFrame(index=[f"{s}_{b}" for b in bar[0]]),
                   var=pd.DataFrame(index=feat[1].astype(str).values))
    a.var_names_make_unique()
    a.obs["sample"] = s; a.obs["group"] = GROUPS[s]
    return a

ads = [load(s) for ss in SAMPLES.values() for s in ss]
A = ad.concat(ads, join="outer", fill_value=0, index_unique=None)
A.var_names_make_unique(); A.obs_names_make_unique()
A.var["mt"] = A.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True)
sc.pp.filter_cells(A, min_genes=200); sc.pp.filter_genes(A, min_cells=3)
A = A[A.obs.n_genes_by_counts < 6000].copy(); A = A[A.obs.pct_counts_mt < 20].copy()
print("QC后:", A.shape)

sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
sc.pp.highly_variable_genes(A, n_top_genes=2000); A = A[:, A.var.highly_variable].copy()

# ---- ComBat 去批次 (按样本) ----
sc.pp.combat(A, key="sample")
print("ComBat 完成")

sc.pp.scale(A, max_value=10); sc.tl.pca(A, n_comps=30)
sc.pp.neighbors(A, n_neighbors=15, n_pcs=30); sc.tl.leiden(A, resolution=1.0, flavor="igraph", n_iterations=2)
print("cluster 数:", A.obs.leiden.nunique())

MARK = {
 'Hepatocyte': ['Alb','Apoa1','Apoa2','Ttr','Serpina1a','Cyp2e1','Cyp2f2'],
 'Cholangiocyte': ['Krt19','Krt7','Sox9','Epcam','Spp1'],
 'Endothelial': ['Pecam1','Cdh5','Stab2','Clec4g','Kdr','Ramp2'],
 'HSC/Mesench': ['Pdgfrb','Des','Acta2','Col1a1','Gfap','Rgs5','Notch3'],
 'Kupffer/Mac': ['Cd68','Clec4f','Vsig4','Marco','Timd4','Csf1r','Adgre1'],
 'Monocyte': ['Ly6c2','Ccr2','Chil3','S100a8','S100a9','Fcnb'],
 'Neutrophil': ['Ly6g','Mpo','Retnlg','S100a8','S100a9'],
 'T cell': ['Cd3e','Cd3d','Cd3g','Cd8a','Cd4','Lck'],
 'NK': ['Nkg7','Klrb1c','Ncr1','Klrd1','Gzma','Prf1'],
 'B cell': ['Cd79a','Cd79b','Ms4a1','Cd19','Ighm'],
 'cDC': ['Itgax','Cd209a','Batf3','Clec9a'],
}
mk = [g for gs in MARK.values() for g in gs if g in A.var_names]
means = A.to_df()[mk].groupby(A.obs.leiden.values).mean()
Z = (means - means.mean()) / (means.std() + 1e-9)
assign = {}
for cl in Z.index:
    best, bs = None, -1e9
    for ct, gs in MARK.items():
        hit = [g for g in gs if g in Z.columns]
        if hit and Z.loc[cl, hit].mean() > bs: bs, best = Z.loc[cl, hit].mean(), ct
    assign[cl] = best
A.obs["celltype"] = A.obs.leiden.map(assign)

ct = A.obs.groupby(["sample", "group", "celltype"], observed=True).size().unstack(fill_value=0)
prop = (ct.div(ct.sum(1), axis=0) * 100)
print("\n每样本细胞比例(%):")
key = ['Kupffer/Mac','Monocyte','HSC/Mesench','Hepatocyte','Neutrophil','T cell','NK','Endothelial','B cell','Cholangiocyte','cDC']
for k in key:
    if k in prop.columns:
        print(f"  {k:13s} " + " ".join(f"{r[1][:4]}:{prop.loc[r,k]:4.1f}" for r in prop.index))
print("\n组均值(%):")
g = A.obs.groupby(["group","celltype"], observed=True).size().unstack(fill_value=0)
print((g.div(g.sum(1),axis=0)*100).round(1).T.to_string())
A.obs[["sample","group","celltype"]].to_csv(f"{OUT}/GSE295293_celltype_combat.csv")
print("\nDONE")
