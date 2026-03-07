# Ulcerative Colitis - MARCO Expression Analysis

Investigating MARCO (macrophage scavenger receptor) expression across single-cell RNA-seq datasets from human and mouse colon tissue.

## Repository Structure

```
├── analysis/                 # Findings, results tables, methodology docs
│   ├── marco_analysis_report.md   # Main report
│   ├── methodology.md             # Data provenance and methods
│   └── *.csv                      # Result tables
├── scripts/
│   └── python/               # Analysis scripts
│       ├── marco_expression.py    # Expression extraction
│       ├── marco_de_pathways.py   # DE + pathway enrichment
│       └── marco_pseudotime.py    # Pseudotime analysis
├── smillie2019/              # Original Smillie et al. 2019 paper code (SCP259)
├── data/                     # Not tracked
│   ├── seurat/               # Seurat .rds objects
│   ├── single-cell-portal/   # SCP259, SCP1845, SCP2038, SCP2760, SCP2771
│   └── gut-atlas/            # Gut Cell Atlas h5ad + metadata
├── objectives/               # Goals, plans, and context docs
└── .gitignore
```

## Datasets

| Source | ID | Species | Tissue | Spatial | Disease Context |
|---|---|---|---|---|---|
| Smillie et al. 2019 | SCP259 | Human | Colon | No | UC |
| Dominguez Conde et al. 2022 | SCP1845 | Human | Multi-organ | No | Normal |
| Avraham-Davidi et al. 2025 | SCP2038 | Mouse | Colon | Slide-seq | Normal/CRC |
| Mayassi et al. 2024 | SCP2760 | Mouse | Colon | No | Normal |
| Mayassi et al. 2024 | SCP2771 | Mouse | Colon | Visium | DSS colitis |
| James et al. 2020 | — | Human | Colon | No | Normal |

## Original Paper

> Intra- and inter-cellular rewiring of the human colon during ulcerative colitis
> Smillie, C.S., Biton, M.B., Ordovas-Montanes J., et al., Cell, 2019.

Original analysis code is preserved in `smillie2019/`.
