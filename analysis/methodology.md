# MARCO Analysis — Methodology & Data Provenance

## Overview

This document describes exactly how each analysis was performed, which datasets were used, and links to the scripts that produced each result.

---

## Scripts

| Script | Purpose | Output Files |
|---|---|---|
| [`scripts/python/marco_expression.py`](../scripts/python/marco_expression.py) | Phase 1-3: Gene detection, expression extraction, cell-type/disease/spatial summaries | `marco_all_datasets_summary.csv`, `scp2771_marco_spatial.csv` |
| [`scripts/python/marco_de_pathways.py`](../scripts/python/marco_de_pathways.py) | DE analysis (MARCO+ vs MARCO-) and GO/KEGG pathway enrichment | `*_marco_de_results.csv`, `*_marco_enrichment.csv` |
| [`scripts/python/marco_pseudotime.py`](../scripts/python/marco_pseudotime.py) | Diffusion pseudotime on myeloid cells | `*_pseudotime.csv` |

---

## Datasets & How They Were Processed

### 1. SCP259 — Smillie et al., *Cell*, 2019
- **Source**: [Single Cell Portal SCP259](https://singlecell.broadinstitute.org/single_cell/study/SCP259)
- **Paper**: DOI: 10.1016/j.cell.2019.06.029
- **Species**: Human | **Tissue**: Colon (lamina propria biopsies)
- **Files used**:
  - Expression: `data/single-cell-portal/SCP259/expression/.../gene_sorted-Imm.matrix.mtx` (genes × cells MTX)
  - Genes: `Imm.genes.tsv` (20,529 genes)
  - Barcodes: `Imm.barcodes2.tsv` (210,614 cells)
  - Metadata: `metadata/all.meta2.txt` — cell type from `Cluster` column, disease from `Health` column
- **Gene searched**: `MARCO` (human)
- **Expression analysis**: Extracted MARCO row from MTX, joined with metadata. Summarized by `Cluster` and `Health`.
- **DE analysis**: Subset to myeloid cell types (Macrophages, Inflammatory Monocytes, Cycling Monocytes, DC1, DC2 — 21,513 cells, 526 MARCO+). Wilcoxon rank-sum test via `scanpy.tl.rank_genes_groups`. Enrichment via `gseapy.enrichr` on significant genes (padj < 0.05, |logFC| > 0.25).
- **Note**: Epithelial compartment (Epi) also checked — 9 MARCO+ cells out of 123,006, negligible. Fibroblast compartment (Fib) has no MARCO gene.

### 2. SCP1845 — Dominguez Conde et al., *Science*, 2022
- **Source**: [Single Cell Portal SCP1845](https://singlecell.broadinstitute.org/single_cell/study/SCP1845)
- **Paper**: DOI: 10.1126/science.abl5197
- **Species**: Human | **Tissue**: 16 organs (including colon, ileum, lung, liver, blood, etc.)
- **Files used**:
  - Expression: `data/single-cell-portal/SCP1845/expression/.../global_normalized_matrix.mtx.gz` (normalized, genes × cells)
  - Genes: `global_features_normalized.tsv.gz`
  - Barcodes: `global_barcodes_normalized.tsv.gz`
  - Metadata: `metadata/global_meta.tsv.gz` — cell type from `Manually_curated_celltype`, organ from `organ__ontology_label`
- **Gene searched**: `MARCO` (human)
- **Expression analysis**: Extracted MARCO row, joined with metadata. Summarized by cell type and organ.
- **DE analysis**: Subset to 11 myeloid cell types (50,082 cells, 25,303 MARCO+). Same Wilcoxon + Enrichr pipeline.
- **Note**: This is a multi-organ dataset. High MARCO+ rates reflect lung alveolar macrophages dominating the myeloid compartment.

### 3. SCP2038 — Avraham-Davidi et al., *eLife*, 2025
- **Source**: [Single Cell Portal SCP2038](https://singlecell.broadinstitute.org/single_cell/study/SCP2038)
- **Paper**: DOI: 10.7554/eLife.104815
- **Species**: Mouse | **Tissue**: Colon (normal controls from CRC study)
- **Files used**:
  - AnnData: `data/single-cell-portal/SCP2038/anndata/scRNAseq.h5ad` (17,512 cells × 31,053 genes)
  - Cell type from `.obs['labels']` column (Epi, Mac, Mono, B, TNK, etc.)
  - Disease state from `.obs['State']` column (all "normal")
- **Gene searched**: `Marco` (mouse)
- **Expression analysis**: Extracted Marco column from h5ad sparse matrix. Only 4 Marco+ cells total (1 Mac, 3 Epi).
- **DE analysis**: Not performed — too few Marco+ cells.

### 4. SCP2760 — Mayassi et al., *Nature*, 2024
- **Source**: [Single Cell Portal SCP2760](https://singlecell.broadinstitute.org/single_cell/study/SCP2760)
- **Paper**: DOI: 10.1038/s41586-024-08216-z
- **Species**: Mouse | **Tissue**: Colon (SPF, GF, FMT conditions)
- **Files used**:
  - Expression: `data/single-cell-portal/SCP2760/expression/.../counts.mtx.gz` (raw counts, genes × cells)
  - Genes: `features.tsv` (gene name in column 2)
  - Barcodes: `barcodes.tsv`
  - Metadata: `metadata/meta.tsv` — cell type from `cell_annotation`, microbiome from `microbe`
- **Gene searched**: `Marco` (mouse)
- **Expression analysis**: 12 Marco+ cells out of 234,490. Summarized by cell type and microbiome condition.
- **DE analysis**: Not performed — too few Marco+ cells.

### 5. SCP2771 — Mayassi et al., *Nature*, 2024
- **Source**: [Single Cell Portal SCP2771](https://singlecell.broadinstitute.org/single_cell/study/SCP2771)
- **Paper**: DOI: 10.1038/s41586-024-08216-z (same paper as SCP2760)
- **Species**: Mouse | **Tissue**: Colon (Visium spatial, DSS colitis model)
- **Files used**:
  - Expression: `data/single-cell-portal/SCP2771/expression/.../matrix.mtx` (raw counts, genes × spots)
  - Genes: `genes.tsv`
  - Barcodes: `barcodes.tsv`
  - Metadata: `metadata/dss_spatial_meta.tsv` — timepoint from `dss_time_point` (D0, D12, D30, D73)
  - Spatial: `cluster/dss_spatial_cluster.tsv` — X/Y coordinates
- **Gene searched**: `Marco` (mouse)
- **Expression analysis**: 320 Marco+ spots, almost all at D12 (315 spots). Summarized by timepoint. Spatial coordinates exported.
- **DE analysis**: Subset to D12 spots only (3,495 spots, 315 Marco+). Wilcoxon rank-sum. Enrichment via Enrichr with mouse gene names uppercased to query human gene set libraries.
- **Note**: Visium spots are ~55μm and capture multiple cells, so this is not single-cell resolution.

### 6. Gut Cell Atlas — James et al., *Nat Immunol*, 2020
- **Source**: [Gut Cell Atlas](https://www.gutcellatlas.org/)
- **Paper**: DOI: 10.1038/s41590-020-0602-z
- **Species**: Human | **Tissue**: Colon (4 regions: Caecum, Transverse, Sigmoid, mLN)
- **Files used**:
  - AnnData: `data/gut-atlas/Colon_cell_atlas.h5ad` (41,650 cells × 18,927 genes, log-normalized)
  - Cell type from `.obs['cell_type']`, region from `.obs['region']`
- **Gene searched**: `MARCO` (human)
- **Expression analysis**: 112 MARCO+ cells. Summarized by cell type and region.
- **DE analysis**: Subset to 7 myeloid types (662 cells, 98 MARCO+). Wilcoxon rank-sum + Enrichr.

---

## Statistical Methods

### Differential Expression
- **Method**: Wilcoxon rank-sum test (Mann-Whitney U) via `scanpy.tl.rank_genes_groups(method='wilcoxon')`
- **Comparison**: MARCO+ (expression > 0) vs MARCO- (expression = 0) within myeloid cells only
- **Significance thresholds**: padj < 0.05 (Benjamini-Hochberg), |logFC| > 0.25
- **Input data**: Normalized expression values (log-normalized for Gut Atlas and SCP1845; raw counts for SCP259 and SCP2771)

### Pathway Enrichment
- **Method**: Over-representation analysis via `gseapy.enrichr`
- **Gene set libraries**: GO Biological Process 2023, GO Molecular Function 2023, GO Cellular Component 2023, KEGG 2021 Human
- **Input**: Significant upregulated and downregulated gene lists (separate analyses)
- **Significance threshold**: Adjusted P-value < 0.05
- **Mouse gene handling**: For SCP2771, mouse gene names were uppercased (e.g., Marco → MARCO) to query human gene set libraries

### Pseudotime Analysis

Diffusion pseudotime (DPT) was used to order myeloid cells along a monocyte-to-macrophage differentiation trajectory and determine where MARCO+ cells fall. Script: [`scripts/python/marco_pseudotime.py`](../scripts/python/marco_pseudotime.py).

#### Cell Selection

Each dataset was subset to myeloid lineage cells only before pseudotime computation:

| Dataset | Myeloid types included | N cells | Cell type column |
|---|---|---|---|
| Gut Atlas | Monocyte, LYVE1 Macrophage, Macrophage, cycling DCs, cDC1, cDC2, pDC | 662 | `cell_type` |
| SCP1845 | Classical monocytes, Nonclassical monocytes, Intermediate macrophages, Intestinal macrophages, Erythrophagocytic macrophages, Alveolar macrophages, Cycling, DC1, DC2, migDC, pDC | 1,081 | `Manually_curated_celltype` |
| SCP259 | Macrophages, Inflammatory Monocytes, Cycling Monocytes, DC1, DC2 | 21,513 | `Cluster` |

SCP1845 was further subset to gut-only organs (jejunal epithelium, lamina propria, ileum, caecum, sigmoid colon, transverse colon, duodenum) to avoid lung alveolar macrophages dominating the trajectory.

#### Preprocessing Pipeline

The following steps were applied sequentially to each myeloid subset:

1. **Normalization check**: The script inspects a 100×100 sample of the expression matrix to determine if values are integer (raw counts) or floating-point (pre-normalized).
   - Gut Atlas: pre-normalized (log-normalized values from h5ad)
   - SCP1845: pre-normalized (downloaded normalized matrix)
   - SCP259: raw counts → normalized to 10,000 counts per cell (`sc.pp.normalize_total(target_sum=1e4)`), then log-transformed (`sc.pp.log1p`)

2. **Highly variable gene selection**: `sc.pp.highly_variable_genes(n_top_genes=2000, subset=True)` — selects the 2,000 most variable genes and subsets the matrix to only these genes. This removes noise from lowly-expressed/invariant genes.

3. **PCA**: `sc.tl.pca(n_comps=30)` — computes 30 principal components on the HVG-subset matrix. The PCA is computed on the centered (but not scaled) expression matrix.

4. **Nearest neighbor graph**: `sc.pp.neighbors(n_neighbors=15, n_pcs=30)` — builds a k-nearest-neighbor graph using Euclidean distance in PCA space. Each cell is connected to its 15 nearest neighbors. This graph is the foundation for UMAP, diffusion maps, and pseudotime.

5. **UMAP**: `sc.tl.umap()` — computes 2D UMAP embedding from the neighbor graph for visualization. Default parameters: `min_dist=0.5`, `spread=1.0`, using the UMAP algorithm by McInnes et al.

6. **Diffusion map**: `sc.tl.diffmap(n_comps=10)` — computes a diffusion map embedding with 10 components. The diffusion map represents the data as a Markov chain on the neighbor graph, where transition probabilities between cells capture gradual transcriptional changes. The eigenvectors of the transition matrix (diffusion components) capture the dominant axes of variation in a way that is robust to noise and preserves global structure — unlike UMAP which prioritizes local structure.

#### Root Cell Selection

DPT requires a root cell to anchor pseudotime at 0. The root was selected as follows:
- From the specified root cell type (monocytes — the presumed starting point of the monocyte-to-macrophage trajectory)
- Within that cell type, the cell with the **lowest MARCO expression** was chosen as root (since MARCO expression may indicate differentiation progress)

Specific root types used:
| Dataset | Root cell type | Root cell barcode |
|---|---|---|
| Gut Atlas | Monocyte | `GAGGTGATCGGCGGTT-1-4861STDY7208411` |
| SCP1845 gut | Classical monocytes | `Pan_T8010319_AACCGCGTCTCGCATC` |
| SCP259 | Inflammatory Monocytes | `N7.LPA.ATTCGGGAGAGGGT` |

#### Diffusion Pseudotime Computation

`sc.tl.dpt(adata)` — computes diffusion pseudotime from the root cell using the diffusion map. The algorithm:
1. Uses the diffusion map components computed above
2. Computes geodesic distances from the root cell through the diffusion space
3. Assigns each cell a pseudotime value (0 = root, increasing = further from root)
4. Cells unreachable from the root receive `inf` pseudotime (excluded from analysis)

DPT was introduced by Haghverdi et al. (2016, *Nature Methods*) and is implemented in scanpy. It is robust to branching trajectories and noise compared to simple ordering methods.

#### MARCO Expression Along Pseudotime

After computing pseudotime, MARCO expression was analyzed along the trajectory:

1. **Binning**: Pseudotime values were divided into 10 equal-frequency bins (deciles) using `pd.qcut`. Each bin contains approximately N/10 cells.

2. **Per-bin statistics**: For each bin, the following were computed:
   - N cells
   - Mean MARCO expression
   - % MARCO+ cells (expression > 0)
   - Mean pseudotime value

3. **MARCO+ vs MARCO- comparison**: Mann-Whitney U test (`scipy.stats.mannwhitneyu`, two-sided) comparing the pseudotime distributions of MARCO+ and MARCO- cells. This tests whether MARCO+ cells are significantly earlier or later in pseudotime.

4. **Cell type ordering**: Mean and median pseudotime computed per cell type to establish the trajectory order (e.g., monocyte → macrophage → tissue-resident).

5. **Disease stratification** (SCP259 only): Among MARCO+ cells, pseudotime was compared across `Health` categories (Healthy, Non-inflamed, Inflamed) to assess whether disease context shifts where MARCO+ cells appear on the trajectory.

#### Output

Each pseudotime analysis produces a CSV file with one row per cell containing:
- `dpt_pseudotime`: diffusion pseudotime value
- `umap_1`, `umap_2`: UMAP coordinates
- `marco_expr`: MARCO expression level
- `marco_status`: MARCO+ or MARCO-
- Cell type label (column name varies by dataset)

#### Limitations

- **DPT assumes a continuous trajectory**: If myeloid differentiation involves discrete jumps or multiple independent origins, DPT may not capture the true topology. PAGA (partition-based graph abstraction) could complement this analysis.
- **Root cell sensitivity**: Pseudotime ordering depends on the choice of root cell. Different root cells within the same monocyte cluster could produce slightly different orderings, though the overall trajectory structure is robust.
- **No RNA velocity**: scVelo-based RNA velocity would provide directional information (confirming monocyte→macrophage rather than the reverse), but requires spliced/unspliced count matrices not available in these downloads. The monocyte→macrophage direction is assumed based on established biology.
- **SCP259 uses raw counts**: Unlike the other datasets which use pre-normalized data, SCP259 was normalized in-script. This means the absolute pseudotime values are not directly comparable across datasets — only the relative ordering and MARCO distribution patterns are meaningful.
- **SCP1845 gut subset is small**: Only 1,081 gut myeloid cells (vs 50,082 total myeloid across all organs). Some rare cell types (e.g., Intermediate macrophages, n=3) may not be reliably positioned on the trajectory.
- **Visium spatial data (SCP2771) excluded**: Pseudotime was not computed on Visium spots because each spot captures ~5-10 cells. Pseudotime requires single-cell resolution to order individual cells along a trajectory.

### Software Versions
- Python 3.9
- scanpy 1.10.3
- gseapy (latest via pip)
- scipy (for mmread, mannwhitneyu)
- pandas, numpy

---

## Output Files

| File | Description |
|---|---|
| `marco_all_datasets_summary.csv` | Expression summary (N cells, mean, % positive, max) across all datasets by grouping variable |
| `marco_analysis_report.md` | Narrative report with all findings |
| `scp2771_marco_spatial.csv` | Spatial coordinates of Marco+ Visium spots |
| `gut_atlas_marco_de_results.csv` | Full DE table: Gut Cell Atlas myeloid MARCO+ vs MARCO- |
| `gut_atlas_marco_enrichment.csv` | GO/KEGG enrichment: Gut Cell Atlas |
| `scp1845_marco_de_results.csv` | Full DE table: SCP1845 myeloid MARCO+ vs MARCO- |
| `scp1845_marco_enrichment.csv` | GO/KEGG enrichment: SCP1845 |
| `scp259_marco_de_results.csv` | Full DE table: SCP259 UC myeloid MARCO+ vs MARCO- |
| `scp259_marco_enrichment.csv` | GO/KEGG enrichment: SCP259 |
| `scp2771_d12_marco_de_results.csv` | Full DE table: SCP2771 D12 Marco+ vs Marco- spots |
| `gut_atlas_pseudotime.csv` | Pseudotime + UMAP + MARCO for Gut Atlas myeloid cells |
| `scp1845_gut_pseudotime.csv` | Pseudotime + UMAP + MARCO for SCP1845 gut myeloid cells |
| `scp259_pseudotime.csv` | Pseudotime + UMAP + MARCO for SCP259 UC myeloid cells |
| `methodology.md` | This file |
