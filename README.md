# Reversal of liver fibrosis: aetiology-specific residual transcriptome — analysis code

Code accompanying a study that asks which parts of the fibrotic transcriptome fail to resolve
after the injury is removed, and whether the answer depends on the aetiology of the injury or on
the way the injury is reversed. The pipeline re-analyses public liver transcriptomic datasets
(toxin-induced fibrosis in mice, diet-induced steatohepatitis in mice, pharmacological reversal,
and paired human biopsies after bariatric surgery), quantifies per-gene resolution, and compares
immune and extracellular-matrix (ECM) gene sets within and across datasets.

## Layout

```
code/
├── scripts/        analysis and plotting scripts (see pipeline below)
├── data/           input data (not distributed; see data/README.md)
├── out/            all generated outputs (created on first run; git-ignored)
├── requirements.txt
├── LICENSE
└── README.md
```

`scripts/_figstyle.py` holds shared plotting style constants and is imported by the plotting scripts.

## Environment

```
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Tested with Python 3.11+ and the version ranges in `requirements.txt`.

## Data

Input files are not distributed with this repository. `data/README.md` lists every required file,
its public accession, and how to obtain or derive it. Place all files in `data/`.

## Running the pipeline

Scripts read from `data/` and write to `out/`. The two phases must be run in order, because the
plotting scripts read intermediate tables produced by the analysis scripts.

### Phase 1 — analysis

| Script | Output in `out/` |
|---|---|
| `01_reversal_compare.py` | signature sizes and reversal fractions per dataset |
| `02_irreversible_core.py` | `reversibility_ratio.csv` |
| `03_consensus_clusters.py` | `consensus_reversibility.csv` |
| `04_replicate_drug.py` | drug-arm reversibility summary |
| `05_kinetics.py` | `kinetics.csv` |
| `06_verify_GSE261310.py` | second toxin model, non-monotonic check |
| `07_specificity.py` | `specificity_GSE199392.csv`, `specificity_GSE261310.csv` |
| `08_third_dataset.py` | third toxin dataset, module residuals |
| `09_GSE135462.py` | GSE135462 recovery and module residuals |
| `10_unbiased_genesets.py` | Hallmark/GO panel results |
| `11_compartment_within_mode.py` | per-mode compartment residuals |
| `12_human_GSE106737.py` | human per-gene-set paired t |
| `13_toxin_drug_GSE241894.py` | toxin-plus-drug arm comparison |
| `14_composition_vs_state.py` | marker deconvolution and abundance correction |
| `15_mcp_counter.py` | MCP-counter abundance residual |
| `16_singlecell_GSE295293.py` | `GSE295293_celltype.csv`, `GSE295293_celltype_prop.csv` |
| `17_singlecell_combat.py` | `GSE295293_celltype_combat.csv` |
| `18_singlecell_percell.py` | `GSE295293_celltype_percell.csv` |

```bash
cd scripts
for s in 0{1..9}_*.py 1{0..8}_*.py; do python "$s"; done
```

### Phase 2 — figures

| Script | Output in `out/` |
|---|---|
| `30_fig1.py` | `figure1.pdf`, `figure1.tif` |
| `31_fig2.py` | `figure2.pdf`, `figure2.tif` |
| `32_fig3.py` | `figure3.pdf`, `figure3.tif` |
| `33_fig4.py` | `figure4.pdf`, `figure4.tif` |
| `41_figsupp.py` | `figureS1.pdf` … `figureS3.pdf` (+ `.tif`) |

```bash
cd scripts
for s in 3{0..3}_*.py 41_*.py; do python "$s"; done
```

Plotting requires the Phase 1 outputs `out/kinetics.csv` (from `05`), `out/reversibility_ratio.csv`
(from `02`), `out/specificity_*.csv` (from `07`), and `out/GSE295293_celltype*.csv` (from `16`–`18`).

## Notes

- Resolution is quantified per gene as `q = |FC(later)| / |FC(injury)|`, computed only for genes in
  each dataset's fibrosis signature (|log2 fold change| >= 1 versus the dataset-matched control).
- Immune and ECM panels are Hallmark and GO Biological Process gene sets from MSigDB 2023.2;
  the exact set identifiers are written by the analysis scripts.
- The single-cell script distinguishes clustering-based, batch-integrated, and per-cell marker
  scoring annotations; because group and batch are nested in this design, integration is reported
  for comparison only.
