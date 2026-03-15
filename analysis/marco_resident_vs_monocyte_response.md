# MARCO+ Cells in Colitis: Tissue-Resident Reprogramming vs Monocyte Infiltration

## The short answer: it's both, but the dominant signal depends on tissue context

---

## Evidence for infiltrating monocyte-derived MARCO+ cells (clearest in UC)

**SCP259 (human UC)** is the strongest case for monocyte origin:
- MARCO+ cells massively upregulate monocyte markers: S100A8 (+15.9 logFC), S100A9 (+20.2), FCN1, VCAN, SELL, PLAC8
- Tissue-resident markers are strongly **down**: C1QA (-9.5 logFC), C1QB (-7.6), C1QC (-6.5), MAF, VSIG4
- On pseudotime, MARCO+ cells in UC sit **late** (PT=0.68 vs 0.35 for MARCO-), consistent with monocytes that have entered tissue and differentiated partway but haven't acquired resident identity
- These are clearly **recently recruited inflammatory monocytes**, not reprogrammed residents

## Evidence for tissue-resident macrophages gaining MARCO during inflammation (clearest in muscularis)

**GSE237862 (mouse muscularis)** tells a more nuanced story:
- In **naive** muscularis, 15.2% of Cx3cr1+ macrophages express Marco — these are constitutive, tissue-resident Marco+ cells
- These Marco+ residents co-express **Timd4** (embryonic/yolk-sac origin marker, p=1.6e-16) and **Mertk** (p=2.9e-9), both hallmarks of long-lived tissue-resident macrophages, not recent monocyte arrivals
- But they also express **Il6, Tnf, Il1b** — pro-inflammatory cytokines you wouldn't expect from a purely anti-inflammatory resident
- During DSS, the Marco+ rate in Cx3cr1+ macrophages actually **drops** from 15.2% to 10.6%, suggesting the resident Marco+ population is being **diluted** by a massive influx of monocyte-derived cells (Cx3cr1+ cells expand 14.5x, from 138 to 2,001)

So in the muscularis, the baseline Marco+ cells **are** tissue-resident (Timd4+/Mertk+), but they already have a pro-inflammatory lean even at steady state — they're not classically anti-inflammatory like the Lyve1+/CD163+ population. During colitis, new monocytes flood in and also express Marco, diluting the resident signal.

## The Gut Atlas (normal human colon) shows the transitional state

- MARCO+ cells co-express **both** resident markers (LYVE1, FOLR2, CD163) **and** monocyte markers (S100A8/A9, FCN1, CD14)
- On pseudotime, MARCO is **bimodal** — high in early monocytes AND in late LYVE1+ tissue-resident macrophages, with a gap in between
- This suggests that in the normal colon, there are two MARCO+ populations coexisting: incoming monocytes and fully differentiated residents

## Summary by context

| Context | Who's expressing MARCO? | Evidence |
|---|---|---|
| **Normal colon** | Tissue-resident macrophages (LYVE1+/FOLR2+) + some incoming monocytes | Bimodal pseudotime, mixed ontogeny markers |
| **Normal muscularis** | Tissue-resident Cx3cr1+ macrophages (Timd4+/Mertk+) | 15.2% Marco+, embryonic origin markers |
| **UC colon** | Predominantly **infiltrating monocytes** | S100A8/A9+++, C1Q---, late pseudotime |
| **DSS muscularis** | Mix — residents persist but diluted by monocyte influx | Marco+ rate drops from 15.2% to 10.6% despite massive Cx3cr1+ expansion |

## The key nuance

The tissue-resident Marco+ macrophages in the muscularis (Timd4+/Mertk+) are **not classically anti-inflammatory to begin with** — they're already Il6-hi/Tnf-hi/Lyve1-lo at baseline, distinct from the anti-inflammatory Lyve1+/CD163+ population. So it's not quite that anti-inflammatory residents are "switching" to pro-inflammatory. Instead:

1. There's a **constitutive pro-inflammatory-leaning resident** Marco+ population (Cx3cr1+/Timd4+) that persists through colitis
2. During colitis, **monocyte-derived cells flood in** and also express Marco — these dominate the signal in mucosal datasets (SCP259) and dilute the resident signal in muscularis (GSE237862)
3. The truly anti-inflammatory residents (Lyve1+/CD163+) **collapse** during colitis (91% reduction) and express very little Marco in either condition

## What the data can't fully resolve

Whether individual resident Marco+ cells are transcriptionally shifting during inflammation, or whether the composition is just changing due to monocyte influx. A lineage tracing experiment (which Stavely et al. did do — they confirmed residents persist) combined with a Marco reporter would be needed to definitively separate those two possibilities.
