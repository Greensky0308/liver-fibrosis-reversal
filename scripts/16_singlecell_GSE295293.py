#!/usr/bin/env python3
"""16_singlecell_GSE295293.py — 单细胞独立复核: TAA 肝纤维化消退中, 免疫细胞退走了吗?
GSE295293: 小鼠肝 FLEX scRNA, control×3 / fibrosis(TAA)×3 / regression×3.
流程: 载入10x -> QC -> 归一 -> 聚类(Leiden) -> marker 注释细胞类型 -> 比较各组细胞比例.
预期(毒物侧): 免疫细胞(库普弗/单核/中性粒)在 fibrosis 高, regression 仍高(没退走).
"""
import os, gzip
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
    M = scipy.io.mmread(f"{SC}/{s}_matrix.mtx.gz").tocsc().T   # cells x genes
    a = ad.AnnData(X=sp.csr_matrix(M), obs=pd.DataFrame(index=[f"{s}_{b}" for b in bar[0]]),
                   var=pd.DataFrame(index=feat[1].astype(str).values))
    a.var_names_make_unique()
    a.obs["sample"] = s; a.obs["group"] = GROUPS[s]
    return a

print("载入样本 ...")
ads = [load(s) for ss in SAMPLES.values() for s in ss]
A = ad.concat(ads, join="outer", fill_value=0, index_unique=None)
A.var_names_make_unique(); A.obs_names_make_unique()
print("总细胞:", A.shape, dict(A.obs.group.value_counts()))

# QC
A.var["mt"] = A.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True)
sc.pp.filter_cells(A, min_genes=200); sc.pp.filter_genes(A, min_cells=3)
A = A[A.obs.n_genes_by_counts < 6000].copy()
A = A[A.obs.pct_counts_mt < 20].copy()
print("QC后:", A.shape, dict(A.obs.group.value_counts()))

sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
sc.pp.highly_variable_genes(A, n_top_genes=2000); A = A[:, A.var.highly_variable].copy()
sc.pp.scale(A, max_value=10); sc.tl.pca(A, n_comps=30)
sc.pp.neighbors(A, n_neighbors=15, n_pcs=30); sc.tl.leiden(A, resolution=1.0, flavor="igraph", n_iterations=2)
print("cluster 数:", A.obs.leiden.nunique())

# marker 注释
MARK = {
 'Hepatocyte': ['Alb','Apoa1','Apoa2','Ttr','Serpina1a','Cyp2e1','Cyp2f2'],
 'Cholangiocyte': ['Krt19','Krt7','Sox9','Epcam','Spp1'],
 'Endothelial': ['Pecam1','Cdh5','Stab2','Clec4g','Kdr','Ramp2'],
 'HSC/Mesench': ['Pdgfrb','Des','Acta2','Col1a1','Gfap','Rgs5','Notch3'],
 'Kupffer/Mac': ['Cd68','Clec4f','Vsig4','Marco','Timd4','Csf1r','Adgre1','Clec4f'],
 'Monocyte': ['Ly6c2','Ccr2','Chil3','S100a8','S100a9','Fcnb'],
 'Neutrophil': ['Ly6g','Mpo','Retnlg','S100a8','S100a9'],
 'T cell': ['Cd3e','Cd3d','Cd3g','Cd8a','Cd4','Lck'],
 'NK': ['Nkg7','Klrb1c','Ncr1','Klrd1','Gzma','Prf1'],
 'B cell': ['Cd79a','Cd79b','Ms4a1','Cd19','Ighm'],
 'cDC': ['Itgax','Cd209a','Batf3','Clec9a'],
}
mk = [g for gs in MARK.values() for g in gs if g in A.var_names]
sc.tl.score_genes(A, mk, score_name="_tmp")   # warm-up
Xc = pd.DataFrame(index=[])
means = A.to_df()[mk].groupby(A.obs.leiden.values).mean()
Z = (means - means.mean()) / (means.std() + 1e-9)
assign = {}
for cl in Z.index:
    best, bs = None, -1e9
    for ct, gs in MARK.items():
        hit = [g for g in gs if g in Z.columns]
        if hit:
            s = Z.loc[cl, hit].mean()
            if s > bs: bs, best = s, ct
    assign[cl] = best
A.obs["celltype"] = A.obs.leiden.map(assign)
print("\ncluster -> celltype:")
for cl, ct in sorted(assign.items(), key=lambda x: int(x[0])):
    print(f"  c{cl:>3s} ({int((A.obs.leiden==cl).sum()):5d} cells) -> {ct}")

# 比例
ct = A.obs.groupby(["group", "celltype"], observed=True).size().unstack(fill_value=0)
prop = ct.div(ct.sum(1), axis=0)
print("\n细胞比例 (%):")
print((prop * 100).round(1).T.to_string())
prop.to_csv(f"{OUT}/GSE295293_celltype_prop.csv")
A.obs[["sample","group","celltype"]].to_csv(f"{OUT}/GSE295293_celltype.csv")
print("\nDONE")
