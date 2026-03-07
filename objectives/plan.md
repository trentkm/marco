## MARCO Expression Analysis Plan (COMPLETED)

### Phase 0: Repo Reorganization
Restructure `/Volumes/repos/uc/` from a flat layout into an organized workspace.

**Moves:**
- R scripts (analysis.r, colors.r, contamination.r, downsample.r, install.r, markers.r, mm_utils.r, mtx.r, parallel.r, plot.r, run.r, scores.r, tpm.r, run_phenograph.py, cell_subsets.txt) → `smillie2019/`
- Seurat objects (train.*.seur.rds) → `data/seurat/`
- `single-cell-data-sets/` → `data/single-cell-portal/`
- `gut_atlas_data/` → `data/gut-atlas/`
- `gut_atlas_handoff.md` → `objectives/`
- Create `analysis/` and `analysis/results/` directories
- Update `.gitignore` for large data files
- Update `README.md` to reflect new structure

**Target layout:**
```
uc/
├── README.md
├── .claude/
├── smillie2019/           # original Smillie et al. 2019 paper code
├── data/
│   ├── seurat/            # .rds files
│   ├── single-cell-portal/ # SCP259, SCP1845, SCP2038, SCP2760, SCP2771
│   └── gut-atlas/         # Gut Cell Atlas h5ad + metadata
├── objectives/            # goals, plans, handoff docs
│   ├── single-cell-goals.md
│   ├── gut_atlas_handoff.md
│   └── plan.md
└── analysis/              # new MARCO analysis
    ├── marco_expression.py
    └── results/
```

### Phase 1: Check for MARCO in Gene Lists
For each dataset, check whether `MARCO` (human) or `Marco` (mouse) appears in the gene/feature files. Quick filter to know which datasets are usable.

**Datasets:**
| ID | Species | Organ | Spatial? | Disease Context |
|---|---|---|---|---|
| SCP259 | Human | Colon | No | UC (inflamed/non-inflamed) |
| SCP1845 | Human | Ileum/Colon | No | Normal |
| SCP2038 | Mouse | Colon | Yes (Slide-seq) | Normal vs disease |
| SCP2760 | Mouse | Colon | No | Normal |
| SCP2771 | Mouse | Colon | Yes (Visium) | DSS colitis |
| Gut Atlas | Human | Colon | No | Normal (already analyzed, see gut_atlas_handoff.md) |

### Phase 2: Extract MARCO Expression
For each dataset where MARCO exists:
- Load expression matrix (MTX + barcodes + genes) or h5ad
- Extract the MARCO row/column
- Join with metadata (cell type, treatment, spatial coords, etc.)

### Phase 3: Analyze & Characterize
- **By cell type**: Which cell types express MARCO?
- **By species**: Compare human vs mouse patterns
- **By disease state**: Expression in UC/DSS vs normal (SCP259, SCP2771, SCP2038)
- **By spatial location**: Map MARCO+ cells in colon (SCP2038, SCP2771)

### Phase 4: Paper Research
- Spawn subagents to look up each SCP study on Single Cell Portal and PubMed
- Document paper context, experimental design, and relevance to MARCO biology

### Prior Work
The Gut Cell Atlas analysis (see `gut_atlas_handoff.md`) already found:
- MARCO is myeloid-restricted (LYVE1 macrophages, monocytes, macrophages)
- Regional variation: caecum monocytes highest (68.4% MARCO+)
- Donor variability: 10-27% of macrophages are MARCO+
- DCs show 0% MARCO expression
