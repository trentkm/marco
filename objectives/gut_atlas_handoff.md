# Gut Cell Atlas - MARCO Analysis Handoff

## Objective
Exploring MARCO (macrophage scavenger receptor) expression in publicly available gut single-cell RNA-seq datasets from the Gut Cell Atlas (https://www.gutcellatlas.org/).

## Downloaded Data

All files are in `gut_atlas_data/`:

| File | Size | Description |
|---|---|---|
| `Colon_cell_atlas.h5ad` | 503 MB | Full Colon Immune Atlas - 41,650 cells x 18,927 genes. Log-normalized counts. |
| `Colon_immune_metadata.csv` | 8.2 MB | Cell metadata with cell_type, donor, region, BCR/TCR info |
| `Colon_immune_gene_names.csv` | 270 KB | Gene names with index positions |

### H5AD Structure
- **adata.obs columns**: donor, region, n_genes, percent_mito, n_counts, cell_type, BCR/TCR chain info
- **Regions**: Caecum, Transverse, Sigmoid, mLN
- **Donors**: 290b, 298c, 302c, 390c, 417c
- **25 cell types** annotated (see findings below)

## Key Findings So Far

### MARCO expression is myeloid-restricted

| Cell Type | n cells | Mean Expr | % Expressing | Max Expr |
|---|---|---|---|---|
| LYVE1 Macrophage | 91 | 0.574 | 35.2% | 3.52 |
| Monocyte | 98 | 0.482 | 36.7% | 2.26 |
| Macrophage | 268 | 0.161 | 10.4% | 3.72 |
| cycling DCs | 47 | 0.048 | 4.3% | 1.95 |
| All lymphocytes/other | ~41K | ~0 | 0% | - |

### Regional variation in MARCO+ myeloid cells

- **Caecum monocytes**: 68.4% MARCO+ (highest)
- **mLN monocytes**: 50.0% MARCO+
- **LYVE1 macrophages**: consistently 35-50% MARCO+ across all regions
- **Conventional macrophages**: 6-20% MARCO+ depending on region
- **DCs (cDC1, cDC2, pDC)**: 0% MARCO expression

### Donor variability
MARCO+ rate among macrophages (Mac + LYVE1 Mac combined) ranges from 10.4% (donor 290b) to 27.8% (donor 390c).

## Additional Datasets Available (Not Yet Downloaded)

### From the Space-Time Gut Cell Atlas (428K cells, fetal/pediatric/adult)
- **Myeloid-specific h5ad** (focused subset, smaller file):
  - Normalized: `https://cellgeni.cog.sanger.ac.uk/gutcellatlas/myeloid_log_counts02_v2.h5ad`
  - Raw: `https://cellgeni.cog.sanger.ac.uk/gutcellatlas/myeloid_raw_counts02_v2.h5ad`
- Also has Epithelium, Mesenchyme, Endothelium, T/NK, B cell subsets available

### From the Pan-GI Cell Atlas (~1.1M cells across entire GI tract)
- Myeloid lineage subset available in H5AD and RDS formats
- All lineages downloadable: `wget https://cellgeni.cog.sanger.ac.uk/gutcellatlas/pangi/objects_pangi_20251029.zip`
- Interactive cellxgene viewers available for each lineage
- Page: https://www.gutcellatlas.org/pangi.html

### Fetal/Pediatric Atlas
- 62,849 fetal cells + 22,500 pediatric Crohn's disease cells
- Could compare MARCO in healthy vs Crohn's context

## Confirmed Genes Present in This Dataset
MARCO (idx 3006), MSR1 (9955), CD163 (13682), MRC1 (12870), CD68 (17782), CSF1R (6547), ITGAM (17183), CD14 (6422), LYVE1 (11683), CD36 (8336)

## Suggested Next Steps
1. **Co-expression analysis**: Check MARCO vs other macrophage markers (CD163, MRC1, LYVE1, CD68) to characterize the MARCO+ population
2. **Download the myeloid-specific h5ad** from the Space-Time atlas for deeper myeloid subtype resolution
3. **Differential expression**: MARCO+ vs MARCO- macrophages to find co-regulated genes
4. **Pan-GI atlas**: Compare MARCO across stomach, small intestine, colon
5. **Crohn's disease dataset**: Check if MARCO expression changes in IBD context
6. **Visualization**: UMAP/dotplots of MARCO and related markers

## Code Notes
- Loaded with `scanpy` (installed via pip3): `sc.read_h5ad('gut_atlas_data/Colon_cell_atlas.h5ad')`
- Matrix is sparse; use `.todense()` when extracting single gene vectors
- The h5ad has old-format neighbor graphs (anndata prints FutureWarning about moving distances/connectivities to .obsp)
