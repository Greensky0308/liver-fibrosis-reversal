#!/usr/bin/env python3
"""18_singlecell_percell.py — 逐细胞 marker 打分直接注释(不聚类、不整合), 作为比例的无偏交叉验证.
理由: 批次(样本)与分组完全嵌套 -> 聚类/整合都可能抹掉组间生物学; 逐细胞打分不受影响.
"""
import os
import numpy as np, pandas as pd
import scipy.io, scipy.sparse as sp
import scanpy as sc, anndata as ad

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SC   = os.path.join(ROOT, "data/GSE295293"); OUT = os.path.join(ROOT, "out")
sc.settings.verbosity = 0
SAMPLES = {
 'control':    ['GSM8944960_Fixed_RNA_control_wt1198','GSM8944961_Fixed_RNA_control_wt1199','GSM8944962_Fixed_RNA_control_wt1200'],
 'fibrosis':   ['GSM8944963_Fixed_RNA_firosis_wt11','GSM8944964_Fixed_RNA_firosis_wt13','GSM8944965_Fixed_RNA_firosis_wt14'],
 'regression': ['GSM8944966_Fixed_RNA_regression_wt1815','GSM8944967_Fixed_RNA_regression_wt1816','GSM8944968_Fixed_RNA_regression_wt1817'],
}
GROUPS = {s: g for g, ss in SAMPLES.items() for s in ss}
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
def load(s):
    feat = pd.read_csv(f"{SC}/{s}_features.tsv.gz", sep="\t", header=None)
    bar = pd.read_csv(f"{SC}/{s}_barcodes.tsv.gz", sep="\t", header=None)
    M = scipy.io.mmread(f"{SC}/{s}_matrix.mtx.gz").tocsc().T
    a = ad.AnnData(X=sp.csr_matrix(M), obs=pd.DataFrame(index=[f"{s}_{b}" for b in bar[0]]),
                   var=pd.DataFrame(index=feat[1].astype(str).values))
    a.var_names_make_unique(); a.obs["sample"] = s; a.obs["group"] = GROUPS[s]
    return a
A = ad.concat([load(s) for ss in SAMPLES.values() for s in ss], join="outer", fill_value=0, index_unique=None)
A.var_names_make_unique(); A.obs_names_make_unique()
A.var["mt"] = A.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True)
sc.pp.filter_cells(A, min_genes=200); sc.pp.filter_genes(A, min_cells=3)
A = A[A.obs.n_genes_by_counts < 6000].copy(); A = A[A.obs.pct_counts_mt < 20].copy()
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
print("QC后细胞:", A.shape)

for ct, gs in MARK.items():
    hit = [g for g in gs if g in A.var_names]
    sc.tl.score_genes(A, hit, score_name=f"sc_{ct}", ctrl_size=50)
scols = [f"sc_{ct}" for ct in MARK]
S = A.obs[scols].copy(); S.columns = list(MARK)
A.obs["celltype"] = S.idxmax(1)
A.obs["maxscore"] = S.max(1)
A.obs["celltype"] = np.where(A.obs["maxscore"] < 0.05, "Unassigned", A.obs["celltype"])
print("未归类比例: %.1f%%" % (100 * (A.obs.celltype == "Unassigned").mean()))

ct = A.obs.groupby(["sample","group","celltype"], observed=True).size().unstack(fill_value=0)
prop = ct.div(ct.sum(1), axis=0) * 100
print("\n每样本细胞比例(%):")
for k in ['Kupffer/Mac','Monocyte','HSC/Mesench','Hepatocyte','Neutrophil','T cell','NK','Endothelial','B cell','Cholangiocyte','cDC']:
    if k in prop.columns:
        print(f"  {k:13s} " + " ".join(f"{r[1][:4]}:{prop.loc[r,k]:4.1f}" for r in prop.index))
g = A.obs.groupby(["group","celltype"], observed=True).size().unstack(fill_value=0)
print("\n组均值(%):"); print((g.div(g.sum(1),axis=0)*100).round(1).T.to_string())
A.obs[["sample","group","celltype","maxscore"]].to_csv(f"{OUT}/GSE295293_celltype_percell.csv")
print("\nDONE")
