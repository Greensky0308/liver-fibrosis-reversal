# Input data

Input files are not distributed with this repository. Download each file from its public source
and place it in this `data/` directory with the exact filename below. All gene-expression
datasets are public on GEO; no individual-level or controlled data are redistributed here.

## Expression datasets (Gene Expression Omnibus)

| Filename | Accession | Description |
|---|---|---|
| `GSE199392_processed_gene_expression_matrix.xlsx` | GSE199392 | Mouse CCl4: control, fibrosis, regression 1/3/6 wk (supplementary processed matrix) |
| `GSE199392_series_matrix.txt.gz` | GSE199392 | Same series, sample metadata (for manual confirmation of group labels; not read by the scripts) |
| `GSE261310_Matrix.GRCm38.xlsx` | GSE261310 | Mouse TAA: control, fibrosis, recovery 1/2/4 wk |
| `GSE135462_m2.xlsx` | GSE135462 | Mouse CCl4 progression/regression supplementary file |
| `GSE135462_FPKM_recovered.csv` | GSE135462 | FPKM matrix recovered from the supplementary file (does not build it) |
| `GSE241894_counts.csv.gz` | GSE241894 | Mouse CCl4 with 4MU / siHAS2 / siLuc arms and controls (counts) |
| `GSE196908_rpkm.txt.gz` | GSE196908 Exp1 | Mouse GAN diet: chow, vehicle, lanifibranor, semaglutide |
| `GSE196908_Exp2_rpkm.txt.gz` | GSE196908 Exp2 | Mouse GAN diet: chow, vehicle, chow-reversal |
| `GSE196908_series_matrix.txt.gz` | GSE196908 | Sample metadata (needed to map columns to groups) |
| `GSE272274_RPKMCodingGenes.xls.gz` | GSE272274 | Mouse high-fat/fructose diet with mARC1-ASO |
| `GSE106737_series_matrix.txt.gz` | GSE106737 | Human paired liver biopsies before/after RYGB |
| `GSE173961_mRNA_FPKM.txt.gz` | GSE173961 | Mouse TAA with recovery |
| `GSE295293/` | GSE295293 | Mouse TAA single-cell (10x) counts |

## Gene-set and annotation resources

| Filename | Source |
|---|---|
| `Mus_musculus.gene_info.gz` | NCBI gene_info for Mus musculus (symbol and Ensembl xrefs) |
| `mh.all.v2023.2.Mm.symbols.gmt` | MSigDB 2023.2 Hallmark, mouse symbols |
| `m5.go.bp.v2023.2.Mm.symbols.gmt` | MSigDB 2023.2 GO Biological Process, mouse symbols |
| `h.all.v2023.2.Hs.symbols.gmt` | MSigDB 2023.2 Hallmark, human symbols |
| `c5.go.bp.v2023.2.Hs.symbols.gmt` | MSigDB 2023.2 GO Biological Process, human symbols |
| `gpl16686_map.tsv` | GEO platform GPL16686 probe-to-symbol annotation (for GSE106737) |

## Prepared helper files

| Filename | How to obtain |
|---|---|
| `mcpc.txt` | MCP-counter gene signatures (population marker lists). Two tab-separated columns: gene symbol and population name. Distributed with the MCP-counter tool. |
| `gse106737_raw.pkl` | A pickle `(header, rows)` parsed from `GSE106737_series_matrix.txt.gz`, where each entry of `rows` is a probe-level record with the PROBEID first followed by one value per sample in `!Sample_title` order. Regenerate it from the series matrix before running the human script. |
