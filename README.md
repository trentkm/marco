# Ulcerative Colitis - MARCO Expression Analysis

Investigating MARCO (macrophage scavenger receptor) expression across single-cell RNA-seq datasets from human and mouse colon tissue.

## Repository Structure

```
├── smillie2019/              # Original Smillie et al. 2019 paper code (SCP259)
├── data/
│   ├── seurat/               # Seurat .rds objects (not tracked)
│   ├── single-cell-portal/   # SCP downloads (not tracked)
│   └── gut-atlas/            # Gut Cell Atlas data (not tracked)
├── objectives/               # Goals, plans, and context docs
├── analysis/                 # MARCO expression analysis scripts
│   └── results/              # Output tables and plots
└── .gitignore
```

## Datasets

| Source | ID | Species | Tissue | Spatial | Disease Context |
|---|---|---|---|---|---|
| Smillie et al. 2019 | SCP259 | Human | Colon | No | UC |
| Elmentaite et al. | SCP1845 | Human | Ileum/Colon | No | Normal |
| Parigi et al. | SCP2038 | Mouse | Colon | Slide-seq | Normal/disease |
| | SCP2760 | Mouse | Colon | No | Normal |
| | SCP2771 | Mouse | Colon | Visium | DSS colitis |
| Gut Cell Atlas | — | Human | Colon | No | Normal |

## Original Paper

> Intra- and inter-cellular rewiring of the human colon during ulcerative colitis
> Smillie, C.S., Biton, M.B., Ordovas-Montanes J., et al., Cell, 2019.

Original analysis code is preserved in `smillie2019/`.
