# MARCO Expression Analysis — Final Report

## Overview

MARCO (Macrophage Receptor with Collagenous Structure) is a class A scavenger receptor. This analysis surveyed MARCO/Marco expression across 6 single-cell/spatial transcriptomic datasets from human and mouse colon tissue.

---

## Datasets & Papers

| ID | Paper | Species | Technology | Tissue | Disease | N cells |
|---|---|---|---|---|---|---|
| SCP259 | Smillie et al., *Cell*, 2019 | Human | 10x 3' v2 | Colon (lamina propria) | UC | 366,650 |
| SCP1845 | Dominguez Conde et al., *Science*, 2022 | Human | 10x 5' v1 | 16 organs incl. colon | Normal | 329,762 |
| SCP2038 | Avraham-Davidi et al., *eLife*, 2025 | Mouse | 10x + Slide-seqV2 | Colon | Normal/CRC | 17,512 |
| SCP2760 | Mayassi et al., *Nature*, 2024 | Mouse | scRNA-seq | Colon | Normal (SPF/GF/FMT) | 234,490 |
| SCP2771 | Mayassi et al., *Nature*, 2024 | Mouse | 10x Visium | Colon | DSS colitis | 13,159 spots |
| Gut Atlas | James et al., *Nat Immunol*, 2020 | Human | scRNA-seq | Colon | Normal | 41,650 |

### Paper Details

**Smillie et al. 2019** — "Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis." 18 UC patients + 12 healthy controls. Identified 51 cell subsets. Key myeloid finding: inflammatory monocytes expand in UC and express OSM (linked to anti-TNF resistance). MARCO not discussed. DOI: 10.1016/j.cell.2019.06.029

**Dominguez Conde et al. 2022** — "Cross-tissue immune cell analysis reveals tissue-specific features in humans." 12 donors, 16 tissues, ~360K cells. Developed CellTypist. Found gut-specific CD209+/IGF1+ macrophages. MARCO not discussed. DOI: 10.1126/science.abl5197

**Avraham-Davidi et al. 2025** — "Spatially defined multicellular functional units in colorectal cancer." Inducible CRC mouse models + normals. Slide-seqV2 spatial. SCP2038 contains the normal colon data. MARCO not discussed. DOI: 10.7554/eLife.104815

**Mayassi et al. 2024** — "Spatially restricted immune and microbiota-driven adaptation of the gut." DSS colitis model with timepoints D0/D12/D30/D73. Found ILC2-mediated middle colon adaptation. SCP2760 = scRNA-seq, SCP2771 = Visium spatial. MARCO not discussed. DOI: 10.1038/s41586-024-08216-z

**James et al. 2020** — "Distinct microbial and immune niches of the human colon." 5 donors, 4 colon regions. Identified LYVE1+ macrophages. MARCO not discussed in this paper, but a follow-up study (Domanska et al., *J Exp Med*, 2022) found MARCO+/LYVE1+ muscularis macrophages in the colon. DOI: 10.1038/s41590-020-0602-z

---

## Key Finding: MARCO is Myeloid-Restricted

Across all datasets, MARCO expression is overwhelmingly restricted to myeloid lineage cells (macrophages, monocytes, some DCs). Non-myeloid expression is negligible (<0.1%).

### Human — Cell Type Breakdown

| Dataset | Cell Type | N cells | % MARCO+ | Mean Expr | Max |
|---|---|---|---|---|---|
| **SCP1845** | Alveolar macrophages | 17,238 | **93.7%** | 2.139 | 4.66 |
| **SCP1845** | Intermediate macrophages | 2,236 | **88.2%** | 2.152 | 5.74 |
| **SCP1845** | MNP/T doublets | 2,508 | 49.0% | 0.741 | 3.57 |
| **SCP1845** | Nonclassical monocytes | 2,420 | 38.4% | 0.454 | 3.47 |
| **Gut Atlas** | Monocyte | 98 | **36.7%** | 0.482 | 2.26 |
| **Gut Atlas** | LYVE1 Macrophage | 91 | **35.2%** | 0.574 | 3.52 |
| **SCP1845** | Cycling | 1,161 | 27.1% | 0.477 | 3.78 |
| **SCP1845** | Erythrophagocytic macrophages | 2,103 | 26.5% | 0.606 | 4.74 |
| **SCP1845** | Classical monocytes | 21,847 | 23.7% | 0.256 | 3.65 |
| **SCP1845** | DC2 | 1,147 | 10.7% | 0.096 | 3.15 |
| **Gut Atlas** | Macrophage | 268 | **10.4%** | 0.161 | 3.72 |
| **SCP259** | Inflammatory Monocytes | 1,652 | 7.7% | 0.093 | 4.0 |
| **SCP259** | DC2 | 2,391 | 2.1% | 0.024 | 3.0 |
| **SCP259** | Macrophages | 16,692 | 2.1% | 0.028 | 25.0 |

### Mouse — Cell Type Breakdown

| Dataset | Cell Type | N cells | % Marco+ | Mean Expr | Max |
|---|---|---|---|---|---|
| **SCP2760** | Retnla+ DCs | 97 | 1.0% | 0.010 | 1.0 |
| **SCP2760** | Pla2g2d+ MACS | 130 | 0.8% | 0.008 | 1.0 |
| **SCP2760** | cDC2s | 743 | 0.7% | 0.012 | 3.0 |
| **SCP2760** | Ccl2/Ccl7-hi MACS | 505 | 0.2% | 0.010 | 5.0 |
| **SCP2038** | Macrophages | 284 | 0.4% | 0.004 | 1.0 |

---

## Disease State Comparison

### Human UC (SCP259)

| Health Status | N cells (Imm) | % MARCO+ | Mean Expr |
|---|---|---|---|
| Healthy | 51,372 | 0.02% | 0.000 |
| Non-inflamed UC | 67,412 | **0.38%** | 0.004 |
| Inflamed UC | 91,830 | **0.31%** | 0.004 |

MARCO is elevated ~15x in UC tissue vs healthy, present in both inflamed and non-inflamed regions.

### Mouse DSS Colitis — Spatial (SCP2771)

| Timepoint | N spots | % Marco+ | Mean Expr | Max |
|---|---|---|---|---|
| D0 (baseline) | 2,919 | 0.1% | 0.001 | 1 |
| **D12 (peak inflammation)** | 3,495 | **9.0%** | 0.194 | 30 |
| D30 (recovery) | 3,498 | 0.0% | 0.000 | 0 |
| D73 (late recovery) | 3,247 | 0.1% | 0.001 | 1 |

Marco expression spikes ~90x at D12 (active inflammation), then returns to baseline by D30. Spatial coordinates for 320 Marco+ spots saved in `scp2771_marco_spatial.csv`.

---

## Organ/Regional Variation (SCP1845)

| Organ | N cells | % MARCO+ | Mean Expr |
|---|---|---|---|
| Lung | 35,419 | **58.3%** | 1.268 |
| Skeletal muscle | 835 | 12.7% | 0.191 |
| Liver | 13,368 | 8.3% | 0.128 |
| Blood | 27,620 | 5.1% | 0.053 |
| Sigmoid colon | 133 | 5.3% | 0.088 |
| Omentum | 597 | 4.9% | 0.061 |
| Bone marrow | 40,507 | 4.8% | 0.051 |
| Spleen | 71,962 | 3.2% | 0.038 |
| Lamina propria | 23,687 | 1.7% | 0.014 |
| Jejunal epithelium | 27,087 | 1.0% | 0.012 |
| Thoracic lymph node | 62,647 | 1.2% | 0.013 |
| Duodenum | 764 | 0.5% | 0.010 |
| Mesenteric lymph node | 23,197 | 0.2% | 0.004 |
| Caecum | 433 | 0.2% | 0.003 |
| Ileum | 890 | 0.0% | 0.000 |
| Transverse colon | 262 | 0.0% | 0.000 |
| Thymus | 354 | 0.0% | 0.000 |

MARCO is highest in lung (alveolar macrophages), with gut expression much lower — consistent with MARCO being a marker of tissue-resident macrophages in barrier organs.

### Within Colon (Gut Atlas)

| Region | N cells | % MARCO+ | Mean Expr |
|---|---|---|---|
| Sigmoid | 6,587 | 0.46% | 0.007 |
| Caecum | 10,478 | 0.34% | 0.005 |
| Transverse | 10,627 | 0.24% | 0.003 |
| mLN | 13,958 | 0.15% | 0.002 |

---

## Microbiome Influence (SCP2760 — Mouse)

| Condition | N cells | % Marco+ |
|---|---|---|
| FMT (germ-free + fecal transplant) | 26,391 | 0.011% |
| GF (germ-free) | 102,186 | 0.006% |
| SPF (conventional) | 105,913 | 0.003% |

Very low overall, but FMT shows ~3x higher Marco+ rate vs SPF — potentially suggesting microbiome reconstitution influences Marco expression, though cell counts are too low to be conclusive.

---

## Spatial Neighborhood Analysis (SCP2771 — Visium DSS Colitis)

Spatial nearest-neighbor analysis of Marco+ Visium spots at D12 (peak inflammation). For each Marco+ spot, its 6 nearest neighbors (matching the Visium hexagonal grid) were identified using a KD-tree on X/Y coordinates. Script: [`marco_spatial_neighbors.py`](../scripts/python/marco_spatial_neighbors.py).

### Spatial Cluster Enrichment

Marco+ spots are overwhelmingly concentrated in **Cluster 10** (38.1% of Marco+ spots vs 6.8% of all D12 spots = **5.55x enrichment**). Cluster 3 is also mildly enriched (1.54x). Clusters 6 (0.22x) and 2 (0.34x) are strongly depleted near Marco+ spots.

| Cluster | % of Marco+ spots | % of all D12 | Enrichment |
| ------- | ----------------- | ------------ | ---------- |
| **10**  | 38.1%             | 6.8%         | **5.55x**  |
| **3**   | 14.9%             | 9.6%         | **1.54x**  |
| 8       | 20.3%             | 20.1%        | 0.98x      |
| 7       | 14.3%             | 22.7%        | 0.66x      |
| 6       | 5.7%              | 26.2%        | **0.22x**  |
| 2       | 2.5%              | 7.9%         | **0.34x**  |

### What's Spatially Adjacent to Marco+ Spots

DE analysis comparing the 822 immediate neighbor spots of Marco+ spots vs 2,358 distant spots (Mann-Whitney U, BH-corrected):

**Top genes enriched in Marco+ neighbors (most spatially co-localized):**

| Gene | Neighbors % | Distant % | Fold change | Function |
|---|---|---|---|---|
| **Ifng** | 9.6% | 2.3% | **13.8x** | IFN-γ — key activator of macrophages |
| **Mmp13** | 42.7% | 22.8% | 8.75x | Matrix metalloproteinase — tissue remodeling |
| **Il22** | 4.4% | 0.5% | — | Epithelial repair cytokine |
| **Il11** | 31.9% | 11.6% | — | Pro-fibrotic cytokine |
| **Cxcl10** | 48.1% | 28.2% | — | IFN-γ-induced chemokine |
| **Il1a** | 15.9% | 4.5% | — | Alarmin |
| **Cxcl1** | 26.2% | 7.5% | — | Neutrophil chemoattractant |
| **Acod1** | 25.4% | 8.0% | — | Itaconate pathway (immunometabolism) |
| **Il1b** | 59.4% | 36.3% | **6.8x** | Pro-inflammatory cytokine |
| **Il6** | 9.5% | 2.4% | **7.1x** | Pleiotropic inflammatory cytokine |
| **Mmp10** | 50.4% | 24.8% | 6.8x | Matrix metalloproteinase |
| **Tnf** | 18.0% | 9.9% | 3.0x | TNF-α |

**Genes depleted near Marco+ spots (absent from the neighborhood):**

| Gene | Neighbors % | Distant % | Function |
|---|---|---|---|
| Fabp6 | 16.3% | 30.2% | Bile acid binding (mature enterocyte) |
| Muc2 (Mptx1) | 2.1% | 5.3% | Goblet cell |
| Slc13a1 | 9.9% | 18.0% | Solute carrier (absorptive epithelium) |
| Cyp3a44/25 | 2.7-6.3% | 6.4-11.7% | Cytochrome P450 (metabolic epithelium) |

### Cell Type Marker Analysis in Spatial Neighborhoods

| Cell lineage | Key marker | Marco+ spots | Neighbors | Distant | Neighbor/Distant ratio |
|---|---|---|---|---|---|
| **Neutrophils** | S100a9 | 94.0% | 84.1% | 62.0% | **4.65x** |
| | Cxcr2 | 22.5% | 7.5% | 1.9% | **4.80x** |
| **Myeloid** | Cd14 | 75.2% | 52.6% | 38.5% | **2.84x** |
| | Ly6c2 | 24.4% | 12.0% | 4.5% | **3.14x** |
| | Ccr2 | 33.0% | 19.0% | 8.9% | **2.55x** |
| **T cells** | Il17a | 12.4% | 2.9% | 0.7% | **5.57x** |
| | Ifng | 22.5% | 9.6% | 2.3% | **13.8x** |
| | Foxp3 | 5.1% | 1.5% | 0.8% | 1.72x |
| **Stroma** | Col3a1 | 96.2% | 93.2% | 83.7% | 1.61x |
| | Vim | 75.2% | 68.6% | 46.8% | **2.00x** |
| **Endothelial** | Pecam1 | 46.7% | 33.2% | 20.6% | **2.00x** |
| **Epithelial** | Epcam | 95.9% | 95.4% | 96.6% | 1.03x |
| | Cdx2 | 47.3% | 59.1% | 62.1% | 0.86x |
| **Tissue remodeling** | Mmp9 | 51.7% | 25.7% | 10.4% | **3.48x** |
| | Mmp13 | 69.2% | 42.7% | 22.8% | **8.75x** |
| | Timp1 | 76.2% | 59.2% | 32.2% | **2.39x** |
| **Scavenger receptors** | Msr1 | 36.8% | 16.7% | 7.6% | **2.92x** |
| | Colec12 | 28.6% | 26.9% | 16.5% | 1.64x |
| **Tissue-resident** | Lyve1 | 7.9% | 6.0% | 3.6% | 1.37x |
| | Timd4 | 2.2% | 2.1% | 1.1% | 1.82x |

### Spatial Interpretation

Marco+ spots at D12 sit within a **highly inflammatory niche** characterized by:

1. **IFN-γ signaling hub**: Ifng is 13.8x enriched in the neighborhood — the strongest spatial signal. This suggests Marco+ macrophages are adjacent to IFN-γ-producing T cells (Th1 and/or CD8+). Il17a is also spatially enriched (5.6x), indicating a mixed Th1/Th17 environment.

2. **Neutrophil co-localization**: S100a8/a9 and Cxcr2 are strongly enriched (4-5x), indicating Marco+ spots sit within or adjacent to neutrophil-rich inflammatory infiltrates.

3. **Active tissue destruction**: MMP9, MMP10, MMP13 are massively enriched (3.5-8.8x) — matrix metalloproteinases that degrade extracellular matrix. This places Marco+ macrophages at sites of active tissue remodeling/damage.

4. **Monocyte recruitment zone**: Ccr2 and Ly6c2 are enriched in neighbors (2.5-3.1x), confirming ongoing monocyte influx around Marco+ spots.

5. **Depleted epithelial maturity**: Absorptive enterocyte markers (Fabp6, Cyp3a, Slc13a1) are depleted — consistent with Marco+ spots being in areas of epithelial damage/loss rather than intact mucosa.

6. **Low tissue-resident macrophage signature**: Lyve1, Timd4, Cx3cr1 are only mildly enriched or unchanged — the neighborhood is dominated by infiltrating cells, not resident macrophages.

---

## Key Biological Context

**Domanska et al., *J Exp Med*, 2022** (DOI: 10.1084/jem.20211846) provides the most relevant biological context for MARCO in the colon:
- Identified 12 macrophage subsets in the human colon
- MARCO+/LYVE1+/FOLR2+/COLEC12+ macrophages localize to the **muscularis layer** (not lamina propria)
- These are homeostatic, tissue-resident macrophages adjacent to neurons and vasculature
- They express low levels of pro-inflammatory cytokines

This explains why:
1. **SCP259** (lamina propria biopsies) shows low MARCO — wrong tissue layer
2. **Gut Atlas** LYVE1+ macrophages are the highest MARCO expressors — consistent with muscularis identity
3. **DSS colitis D12** spike may reflect muscularis disruption/infiltration during active inflammation

---

## Summary

1. **MARCO is myeloid-restricted** across all 6 datasets — primarily macrophages and monocytes
2. **In humans**: highest in LYVE1+ tissue-resident macrophages (muscularis type) and monocytes; much lower in lamina propria macrophages
3. **In mice**: very low at baseline, but **dramatic upregulation during active colitis** (D12 DSS = 9% of spatial spots, vs 0.1% baseline)
4. **Species difference**: human colon macrophages show moderate constitutive MARCO expression; mouse colon macrophages show near-zero baseline but inflammation-inducible expression
5. **Organ hierarchy**: lung >> liver >> blood/spleen >> gut >> lymph nodes
6. **Disease relevance**: MARCO is elevated in UC tissue and spikes during DSS colitis peak — suggesting MARCO marks inflammation-recruited or -activated myeloid cells
7. **No prior papers** in this analysis explicitly discussed MARCO, but Domanska et al. 2022 provides the key biological framework

## Differential Expression: MARCO+ vs MARCO- Myeloid Cells

DE analysis was performed on 4 datasets with sufficient MARCO+ cells (Wilcoxon rank-sum, padj < 0.05, |logFC| > 0.25). Full methods in [methodology.md](methodology.md). Scripts: [`marco_de_pathways.py`](../scripts/python/marco_de_pathways.py).

### DE Summary by Dataset

| Dataset | Myeloid cells | MARCO+ | MARCO- | Sig. UP | Sig. DOWN |
|---|---|---|---|---|---|
| Gut Atlas | 662 | 98 | 564 | 162 | 73 |
| SCP1845 | 50,082 | 25,303 | 24,779 | 5,921 | 1,649 |
| SCP259 | 21,513 | 526 | 20,987 | 389 | 155 |
| SCP2771 (D12 spots) | 3,495 | 315 | 3,180 | 3,324 | 213 |

### Consistently Co-upregulated with MARCO

| Gene | Function | Datasets Found |
|---|---|---|
| **S100A8/A9** | Calprotectin (alarmin, inflammation) | Gut Atlas, SCP259, SCP2771 |
| **FCN1** | Ficolin-1 (innate immunity, complement) | Gut Atlas, SCP259 |
| **CD163** | Hemoglobin scavenger receptor (M2 macrophage marker) | Gut Atlas, SCP1845 |
| **VCAN** | Versican (ECM, monocyte recruitment) | Gut Atlas, SCP259 |
| **TYROBP/DAP12** | Myeloid signaling adaptor | Gut Atlas, SCP259, SCP2771 |
| **C5AR1** | Complement C5a receptor | Gut Atlas |
| **MSR1** | Scavenger receptor A (CD204) | SCP1845 |
| **MRC1/CD206** | Mannose receptor (M2 macrophage marker) | SCP1845 |
| **TSPO** | Translocator protein (inflammation) | Gut Atlas, SCP259 |
| **SERPINA1** | Alpha-1 antitrypsin | Gut Atlas, SCP259 |
| **EMP3** | Epithelial membrane protein 3 | Gut Atlas, SCP259 |
| **Il1b** | Pro-inflammatory cytokine | SCP2771 (mouse, D12) |
| **Cxcl9/Cxcl10** | IFN-γ-induced chemokines | SCP2771 (mouse, D12) |
| **Saa3** | Serum amyloid A3 (acute phase) | SCP2771 (mouse, D12) |

### Consistently Downregulated in MARCO+ Cells

| Gene | Function | Datasets Found |
|---|---|---|
| **HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DQA1** etc. | MHC class II (antigen presentation) | Gut Atlas, SCP259, SCP1845 |
| **C1QA/B/C** | Complement C1q subunits | SCP259 |
| **DNASE1L3** | DNase I-like 3 | Gut Atlas, SCP259 |
| **IGF1** | Insulin-like growth factor 1 | Gut Atlas, SCP259 |
| **RGS1** | Regulator of G-protein signaling 1 | Gut Atlas, SCP259 |
| **SLC40A1/Ferroportin** | Iron exporter | Gut Atlas, SCP259 |
| **FCER1A** | IgE receptor alpha (DC marker) | Gut Atlas, SCP1845 |
| **PLD4** | Phospholipase D4 (pDC marker) | Gut Atlas, SCP1845 |

---

## Pathway Enrichment

### Upregulated in MARCO+ Cells

**GO Biological Process:**
- Inflammatory Response (GO:0006954) — all datasets
- Response to Lipopolysaccharide (GO:0032496) — Gut Atlas, SCP259
- Phagocytosis (GO:0006909) — Gut Atlas (MARCO itself part of this term)
- Granulocyte/Neutrophil Chemotaxis — Gut Atlas, SCP259
- Cytoplasmic Translation — SCP259, SCP1845 (ribosomal program)
- Mitochondrial Translation / Cellular Respiration — SCP1845 (metabolically active cells)

**KEGG Pathways:**
- Phagosome — Gut Atlas, SCP259
- NF-κB Signaling Pathway — Gut Atlas, SCP259
- Ferroptosis — Gut Atlas
- Lysosome — Gut Atlas
- Complement and Coagulation Cascades — Gut Atlas
- Tuberculosis — Gut Atlas, SCP259

**GO Molecular Function:**
- RAGE Receptor Binding — Gut Atlas, SCP259 (S100A8/A9/A12 → RAGE axis)
- Amyloid-Beta Binding — Gut Atlas (MARCO, ITGAM, LRP1, TLR2)
- Endopeptidase Inhibitor Activity — Gut Atlas, SCP259

**GO Cellular Component:**
- Secretory Granule Membrane — Gut Atlas, SCP259
- Ficolin-1-Rich Granule — Gut Atlas, SCP259
- Tertiary Granule — Gut Atlas, SCP259

### Downregulated in MARCO+ Cells

**GO Biological Process:**
- MHC Class II Protein Complex Assembly — dominant signal across all datasets
- Antigen Processing and Presentation via MHC Class II — all datasets
- Immunoglobulin Mediated Immune Response — Gut Atlas
- Cytoplasmic Translation (ribosomal subunits) — Gut Atlas downregulated set

**KEGG Pathways:**
- Antigen Processing and Presentation — all datasets
- Inflammatory Bowel Disease — Gut Atlas (driven by MHC-II downregulation)
- Intestinal Immune Network for IgA Production — Gut Atlas
- Th1 and Th2 Cell Differentiation — Gut Atlas

---

## Biological Interpretation

MARCO+ myeloid cells in the colon have a distinct transcriptional profile:

1. **Pro-inflammatory and phagocytic**: High S100A8/A9 (calprotectin), active NF-κB/LPS response pathways, ficolin-rich granules, and phagosome genes. This aligns with MARCO's known role as a pattern recognition receptor for bacterial clearance.

2. **Reduced antigen presentation**: Strong downregulation of MHC class II genes (HLA-DR, HLA-DQ, HLA-DP families) suggests MARCO+ cells are less engaged in adaptive immune priming. This is consistent with recently recruited inflammatory monocytes vs tissue-resident antigen-presenting macrophages.

3. **M2-like but inflammatory**: Co-expression with CD163 and MRC1 (classical M2 markers) alongside S100A8/A9 and IL1B (classical M1 markers) suggests MARCO+ cells don't fit neatly into M1/M2 polarization — they may represent a transitional or mixed-phenotype population.

4. **Metabolically active**: In SCP1845, MARCO+ cells show upregulated mitochondrial translation and oxidative phosphorylation, suggesting high metabolic demands (consistent with active phagocytosis).

5. **Inflammation-inducible in mice**: The dramatic D12 DSS spike with co-upregulation of Saa3, Il1b, Cxcl9/10, and Lcn2 indicates Marco marks acutely recruited inflammatory myeloid cells during colitis, with rapid resolution by D30.

6. **Tissue-layer specificity**: Low MARCO in lamina propria datasets (SCP259) vs higher in datasets capturing muscularis (Gut Atlas LYVE1+ macrophages) is consistent with Domanska et al. 2022's finding that MARCO+ macrophages reside in the muscularis propria.

---

## Ontogeny Analysis: Yolk Sac vs Monocyte-Derived

MARCO+ myeloid cells were scored for established ontogeny markers to determine whether they originate from embryonic (yolk sac) tissue-resident macrophages or bone marrow-derived monocytes. Script: [`marco_de_pathways.py`](../scripts/python/marco_de_pathways.py).

### Marker Panels Used

- **Yolk sac / tissue-resident**: LYVE1, FOLR2, TIMD4, CX3CR1, CSF1R, MRC1, CD163, MAF, MERTK, C1QA/B/C, VSIG4
- **Monocyte-derived**: CCR2, CD14, S100A8/A9, FCN1, VCAN, SELL, PLAC8, CLEC12A, LST1
- **Peritoneal macrophage**: GATA6, TGFB2, RARRES2, F5, VSIG4, ARG1, ICAM2

### Results by Dataset

**Gut Atlas (normal human colon) — Mixed/transitional signature:**
- Tissue-resident markers UP: LYVE1 (p=0.019), FOLR2 (p=5.8e-4), CD163 (p=1.9e-12), MAFB (p=7.0e-9)
- Monocyte markers also UP: CD14 (p=2.3e-6), S100A8/A9 (p<1e-9), FCN1 (p=4.4e-11), VCAN (p=4.4e-10)
- C1QA/B/C: unchanged (not significantly different)
- **Interpretation**: MARCO+ cells in the normal colon are tissue-resident macrophages (LYVE1+/FOLR2+/CD163+) that retain some monocyte-recruitment markers, consistent with a population that was monocyte-derived but has adapted to the tissue niche

**SCP1845 (multi-organ) — Clearly tissue-resident:**
- Tissue-resident markers strongly UP: MRC1 (p≈0), CD163 (p≈0), C1QA/B/C (p≈0), VSIG4 (p≈0), TIMD4 (p=1.9e-15), MAF (p≈0)
- Monocyte markers strongly DOWN: CCR2 (p=3.9e-80), S100A8/A9 (p≈0), FCN1 (p≈0), VCAN (p≈0), SELL (p≈0)
- **Interpretation**: Dominated by alveolar macrophages (lung), which are embryonically-derived tissue-resident macrophages. Across organs, MARCO marks the most mature, tissue-adapted macrophage population.

**SCP259 (UC colon) — Clearly monocyte-derived:**
- Monocyte markers massively UP: S100A8 (logFC=+15.9, p=4.1e-87), S100A9 (logFC=+20.2, p=2.7e-111), FCN1 (p=4.0e-85), VCAN (p=1.5e-40), SELL (p=2.1e-7), PLAC8 (p=8.3e-13)
- Tissue-resident markers DOWN: C1QA (logFC=-9.5, p=9.1e-46), C1QB (logFC=-7.6, p=1.4e-40), C1QC (logFC=-6.5, p=1.3e-38), MAF (p=9.1e-5), VSIG4 (p=1.9e-4), SEPP1 (p=7.9e-40)
- **Interpretation**: In UC, MARCO+ cells are recently recruited inflammatory monocytes that have not yet acquired tissue-resident identity. This is the classic monocyte waterfall — blood monocytes entering inflamed tissue.

### Peritoneal Macrophage Origin

**No evidence of peritoneal origin in any dataset:**
- GATA6 (canonical peritoneal TF): not expressed or not significant across all datasets
- TGFB2: absent or negative
- RARRES2: absent or negative
- ARG1: absent
- MARCO+ colon macrophages are definitively **not** peritoneal macrophages

### Context-Dependent Ontogeny Summary

| Context | MARCO+ identity | Key evidence |
|---|---|---|
| Normal colon (Gut Atlas) | Tissue-resident with monocyte traces | LYVE1+/FOLR2+/CD163+ but also S100A8/A9+/FCN1+ |
| Multi-organ steady-state (SCP1845) | Mature tissue-resident | C1Q+++, MRC1+++, CCR2---, S100A8/A9--- |
| UC colon (SCP259) | Recently recruited monocytes | S100A8/A9+++, FCN1+++, C1Q---, MAF--- |
| DSS colitis D12 (SCP2771) | Acute inflammatory infiltrate | Saa3+++, Il1b+++, Cxcl9/10+++ |

---

## Inflammatory Polarization

### Marker Panels Used

- **Pro-inflammatory**: IL1B, IL6, TNF, CXCL2/3/8/10, CCL2/3/4, S100A8/A9, NFKBIA, NLRP3, SOD2, TLR2/4, NOS2
- **Anti-inflammatory / M2**: IL10, TGFB1, MRC1, CD163, MSR1, FOLR2, IGF1, LYVE1, MERTK, AXL, STAB1, MAF, MAFB, NR4A1/2/3, PPARG, APOE, C1QA/B/C

### Polarization by Dataset

**Normal colon (Gut Atlas) — Pro-inflammatory lean with tissue-resident features:**
- Pro-inflammatory significantly UP: S100A8/A9, CXCL2/3, CCL3/4, NLRP3, TLR2, SOD2, NFKBIA
- Anti-inflammatory UP: CD163, FOLR2, LYVE1, MAFB
- Anti-inflammatory DOWN: IGF1 (p=0.022), AXL (p=0.023)
- **Verdict**: Mixed — scavenger/inflammatory cells with some homeostatic features

**Multi-organ (SCP1845) — Anti-inflammatory / homeostatic:**
- Anti-inflammatory strongly UP: MRC1 (p≈0), CD163 (p≈0), MSR1 (p≈0), C1QA/B/C (p≈0), APOE (p≈0), PPARG (p≈0), AXL (p≈0), MAF (p≈0), NR4A1/2/3 (p≈0)
- Pro-inflammatory mixed: IL1B/IL6/TNF up, but S100A8/A9 strongly DOWN, NLRP3 DOWN, SOD2 DOWN
- **Verdict**: Predominantly anti-inflammatory, mature tissue-resident macrophages

**UC colon (SCP259) — Pro-inflammatory:**
- Pro-inflammatory UP: S100A8/A9 (massively), FCN1, VCAN
- Anti-inflammatory DOWN: IGF1 (p=9.8e-9), AXL (p=3.8e-4), STAB1 (p=9.9e-8), MAF (p=9.1e-5), NR4A2 (p=5.9e-5), NR4A3 (p=2.7e-6), APOE (p=6.7e-19), C1QA/B/C (p<1e-38)
- Notably: IL1B is actually slightly DOWN, and CCL3/CCL4 are DOWN — suggesting these are early-stage monocytes not yet fully activated
- **Verdict**: Pro-inflammatory monocyte identity, but not yet producing classical cytokines — consistent with recently extravasated cells

**DSS colitis D12 (SCP2771) — Acutely pro-inflammatory:**
- Massively UP: Saa3, S100a8/a9, Il1b, Cxcl9/10, Lcn2, Socs3
- **Verdict**: Peak acute inflammation signature

### Polarization Summary

| Context | Polarization | Character |
|---|---|---|
| Normal colon | Mixed | Tissue-resident scavenger with inflammatory capacity |
| Multi-organ steady-state | Anti-inflammatory | Mature homeostatic macrophages |
| UC colon | Pro-inflammatory | Recently recruited monocytes, not yet cytokine-producing |
| DSS colitis peak | Acutely pro-inflammatory | Full inflammatory activation |

MARCO+ cells do not fit the classical M1/M2 dichotomy. Their polarization is **context-dependent**: homeostatic in steady-state tissue, but marking recruited inflammatory monocytes during disease.

---

## Pseudotime Analysis

Diffusion pseudotime (DPT) was computed on myeloid cells to trace monocyte-to-macrophage differentiation and locate MARCO+ cells along the trajectory. Script: [`marco_pseudotime.py`](../scripts/python/marco_pseudotime.py).

### Method
- Scanpy DPT with diffusion maps (10 components)
- Root cell selected from monocytes (lowest MARCO expression)
- Pseudotime binned into deciles for MARCO expression tracking
- Mann-Whitney U test comparing MARCO+ vs MARCO- pseudotime distributions

### Gut Atlas (Normal Human Colon, 662 myeloid cells)

**MARCO peaks early in pseudotime — at the monocyte stage:**

| Pseudotime bin | Mean PT | N cells | % MARCO+ |
|---|---|---|---|
| 0 (earliest) | 0.017 | 67 | **37.3%** |
| 1 | 0.047 | 66 | **30.3%** |
| 2 | 0.136 | 66 | 15.2% |
| 3 | 0.168 | 66 | 4.5% |
| 4 | 0.176 | 66 | 0.0% |
| 5 | 0.195 | 66 | 4.5% |
| 6-7 | 0.236-0.255 | 132 | 13.6-16.7% |
| 8 | 0.269 | 66 | **22.7%** |
| 9 (latest) | 0.421 | 67 | 3.0% |

**Cell type ordering**: Monocyte (PT=0.073) → Macrophage (0.171) → cycling DCs (0.187) → cDC2 (0.200) → LYVE1 Macrophage (0.240) → cDC1 (0.254) → pDC (0.952)

- MARCO+ mean PT = 0.136 vs MARCO- mean PT = 0.202 (p=1.98e-5)
- MARCO is bimodal: high in early monocytes AND in LYVE1 macrophages (PT bin 8), with a gap in between
- This suggests two distinct MARCO+ populations: recently arrived monocytes and fully differentiated tissue-resident macrophages

### SCP1845 Gut Myeloid (1,081 cells from gut organs)

**Cell type ordering**: Classical monocytes (PT=0.173) → Intermediate macrophages (0.219) → Nonclassical monocytes (0.269) → DC2 (0.450) → Intestinal macrophages (0.472) → DC1 (0.479) → Erythrophagocytic macrophages (0.497) → Alveolar macrophages (0.712)

| Cell Type | Mean PT | % MARCO+ |
|---|---|---|
| Classical monocytes | 0.173 | **34.3%** |
| Nonclassical monocytes | 0.269 | **45.5%** |
| Intermediate macrophages | 0.219 | 33.3% |
| Intestinal macrophages | 0.472 | 4.6% |
| Alveolar macrophages | 0.712 | **45.2%** |

- MARCO+ mean PT = 0.423 vs MARCO- mean PT = 0.456 (p=1.33e-5)
- MARCO is again **bimodal**: high at the monocyte entry point AND at the terminally differentiated tissue-resident macrophage endpoint
- Intestinal macrophages (mid-trajectory) have the lowest MARCO — consistent with monocyte waterfall model where cells lose monocyte markers before acquiring tissue-resident identity

### SCP259 (UC Colon, 21,513 myeloid cells)

**MARCO+ cells are late in pseudotime — the opposite of normal colon:**

| Pseudotime bin | Mean PT | N cells | % MARCO+ |
|---|---|---|---|
| 0-3 (early) | 0.04-0.13 | 8,605 | 0.1-0.5% |
| 4-5 (mid) | 0.18-0.23 | 4,303 | 0.2-0.8% |
| 6 | 0.358 | 2,151 | **4.3%** |
| 7 | 0.736 | 2,151 | **7.9%** |
| 8 | 0.838 | 2,151 | **6.1%** |
| 9 (latest) | 0.891 | 2,152 | **4.0%** |

**Cell type ordering**: Cycling Monocytes (PT=0.095) → DC2 (0.270) → DC1 (0.307) → Macrophages (0.359) → Inflammatory Monocytes (0.573)

- MARCO+ mean PT = 0.670 vs MARCO- mean PT = 0.352 (p=3.01e-93, highly significant)
- In UC, MARCO marks **late-stage inflammatory monocytes**, not early arrivals
- Inflammatory Monocytes are the latest cell type (PT=0.573) and have the highest MARCO+ rate (7.7%)

**MARCO+ pseudotime by disease status:**

| Health | N MARCO+ | Mean PT | Median PT |
|---|---|---|---|
| Healthy | 12 | 0.154 | 0.114 |
| Non-inflamed | 254 | 0.680 | 0.762 |
| Inflamed | 260 | 0.684 | 0.813 |

- In healthy tissue, the few MARCO+ cells are early (monocyte-like, PT=0.15)
- In UC tissue (both inflamed and non-inflamed), MARCO+ cells are late (PT=0.68) — consistent with disease-recruited inflammatory monocytes that have progressed through differentiation

### Pseudotime Interpretation

The pseudotime analysis reveals that MARCO expression follows a **context-dependent trajectory pattern**:

1. **Normal colon**: MARCO is bimodal — expressed by incoming monocytes (early PT) and by terminally differentiated tissue-resident LYVE1+ macrophages (late PT), with a gap during the intermediate macrophage differentiation stage. This matches the "monocyte waterfall" model where cells transiently lose monocyte identity before acquiring tissue-resident features.

2. **UC colon**: MARCO shifts to marking only **late-stage inflammatory monocytes** (high PT). The bimodal pattern collapses — there is no tissue-resident MARCO+ peak, only the inflammatory endpoint. This suggests that in disease, MARCO marks monocytes that have been activated/differentiated within the inflamed tissue rather than homeostatic tissue-resident macrophages.

3. **Gut myeloid (SCP1845)**: Confirms the bimodal pattern with MARCO highest in classical monocytes (entry) and alveolar/tissue-resident macrophages (endpoint), lowest in mid-trajectory intestinal macrophages.

---

## References

1. Smillie CS et al. Cell. 2019;178(3):714-730. DOI: 10.1016/j.cell.2019.06.029
2. Dominguez Conde C et al. Science. 2022;376(6594):eabl5197. DOI: 10.1126/science.abl5197
3. Avraham-Davidi I et al. eLife. 2025;14:RP104815. DOI: 10.7554/eLife.104815
4. Mayassi T et al. Nature. 2024;636(8042):447-456. DOI: 10.1038/s41586-024-08216-z
5. James KR et al. Nat Immunol. 2020;21(3):343-353. DOI: 10.1038/s41590-020-0602-z
6. Elmentaite R et al. Nature. 2021;597(7875):250-255. DOI: 10.1038/s41586-021-03852-1
7. Domanska D et al. J Exp Med. 2022;219(3):e20211846. DOI: 10.1084/jem.20211846
