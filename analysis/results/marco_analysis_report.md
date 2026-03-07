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

DE analysis was performed on 4 datasets with sufficient MARCO+ cells (Wilcoxon rank-sum, padj < 0.05, |logFC| > 0.25). Full methods in [methodology.md](methodology.md). Scripts: [`marco_de_pathways.py`](../../analysis/marco_de_pathways.py).

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

## References

1. Smillie CS et al. Cell. 2019;178(3):714-730. DOI: 10.1016/j.cell.2019.06.029
2. Dominguez Conde C et al. Science. 2022;376(6594):eabl5197. DOI: 10.1126/science.abl5197
3. Avraham-Davidi I et al. eLife. 2025;14:RP104815. DOI: 10.7554/eLife.104815
4. Mayassi T et al. Nature. 2024;636(8042):447-456. DOI: 10.1038/s41586-024-08216-z
5. James KR et al. Nat Immunol. 2020;21(3):343-353. DOI: 10.1038/s41590-020-0602-z
6. Elmentaite R et al. Nature. 2021;597(7875):250-255. DOI: 10.1038/s41586-021-03852-1
7. Domanska D et al. J Exp Med. 2022;219(3):e20211846. DOI: 10.1084/jem.20211846
