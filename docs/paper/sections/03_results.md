# Results

## AD genetic risk maps to myeloid-rich regions across the whole body

We projected AD GWAS signals onto whole-body spatial transcriptomics data from 6-week-old C57BL/6J mice using gsMap. The resulting spatial risk map revealed significant AD heritability enrichment across all 16 profiled organs, with the strongest signals in immune-rich regions (**Figure 1A–B**). Organ-level Cauchy combination analysis identified Spleen (p = 2.5 × 10⁻⁸), Lymph Node (p = 5.5 × 10⁻⁷), and Lung (p = 5.1 × 10⁻⁶) as the most enriched organs, while Brain ranked lower (p = 4.2 × 10⁻³ to 2.4 × 10⁻²) (**Figure 1C**). Results were highly consistent between biological replicates (Pearson r = 0.94 for annotation-level Cauchy p-values).

## Tissue-resident macrophages universally carry AD risk

Cell-type-level analysis revealed that myeloid/macrophage annotations were the most significantly enriched cell type in 11 of 15 organs examined (**Figure 2A**). Tissue-resident macrophages showed convergent AD enrichment regardless of their organ of residence: Kupffer cells in liver (p = 4.9 × 10⁻¹⁰), alveolar macrophages in lung (p = 1.2 × 10⁻¹⁰), and monocytes in bone marrow (p = 5.2 × 10⁻¹⁰). In brain, microglia — the brain-resident macrophage — was the only significantly enriched cell type (p = 2.4 × 10⁻⁸ after microglia enrichment), while all neuronal subtypes were non-significant (p > 0.85) (**Figure 2B**). This pattern indicates that AD genetic risk operates through the myeloid cell lineage rather than through organ-specific mechanisms.

## Colonic epithelial cells harbor an independent AD risk program

In contrast to the myeloid-dominated pattern in most organs, three gastrointestinal organs — Colon, Ileum, and Stomach — exhibited significant enrichment in non-immune cell types (**Figure 3A**). In Colon, secretory cells (p = 7.3 × 10⁻⁷) and enterocytes (p = 3.8 × 10⁻⁶) ranked above immune cells in enrichment strength. Gene-level analysis revealed that this epithelial enrichment was driven by a distinct set of AD risk genes centered on the amyloid precursor protein processing pathway: APP (PCC = 0.70), SORL1 (PCC = 0.72), PSEN1 (PCC = 0.69), ADAM10 (PCC = 0.71), CD2AP (PCC = 0.82), and PICALM (PCC = 0.73) all showed their highest or near-highest cross-organ correlation with AD risk in Colon (**Figure 3B**). These APP pathway genes were specifically expressed in secretory cells and enterocytes, not in colonic immune cells (**Figure 3C**). Spatial expression mapping confirmed that App, Sorl1, Psen1, and Cd2ap were broadly expressed across colon epithelial regions (**Figure 3D**).

Colon was the only organ where the APP processing gene set exceeded the myeloid gene set in mean AD correlation (0.68 vs 0.27), representing a qualitative reversal of the pattern seen in all other organs (**Figure 3E**). This suggests that colonic epithelium operates an intrinsic APP processing program that is genetically linked to AD risk independently of the immune/myeloid pathway.

## AD risk distribution is stable across the mouse lifespan

To examine whether AD genetic risk changes with aging, we mapped single-cell RNA-seq data from eight sources (covering 15 organs, 10 age points from 1 to 30 months) onto the spatial reference using CMAP. gsMap analysis of these 126 age-resolved pseudo-spatial samples confirmed that the organ enrichment pattern was stable across the lifespan (**Figure 4A**). Spearman correlation analysis showed no significant age trend for any organ (all ρ < 0.5, p > 0.05) (**Figure 4B**). The whole-body spatial enrichment pattern at 3 months was visually indistinguishable from that at 21 months (**Figure 4C**). This age-independence is consistent with AD genetic risk being encoded in the germline rather than acquired through age-dependent epigenetic or transcriptional changes.

## Multi-trait comparison reveals AD-specific peripheral enrichment

To determine whether the observed patterns are AD-specific, we applied the same gsMap pipeline to GWAS data for PD, ALS, MS, and FTD. AD showed dramatically stronger enrichment than all other traits across every organ, with -log₁₀(p) values of 7–10 compared to 3–6 for ALS and 1–4 for PD (**Figure 5A**). 

Myeloid cell enrichment followed a clear disease-specific gradient: AD ≫ ALS > PD (**Figure 5B**). In brain, AD and ALS enrichment was microglia-driven, while PD enrichment was predominantly in oligodendrocytes and astrocytes (p = 6.0 for oligodendrocyte vs p = 3.0 for microglia), representing a fundamentally different cellular mechanism (**Figure 5C**).

The gut epithelial enrichment pattern was the most disease-discriminating feature. In Colon, AD showed strong epithelial enrichment (-log₁₀p = 6.2 for non-immune cells), ALS showed moderate enrichment (-log₁₀p = 4.5 for secretory cells), while PD showed negligible enrichment (-log₁₀p = 2.6) (**Figure 5D**). This indicates that gut epithelial involvement in neurodegenerative disease genetic risk is not a generic feature but is preferentially associated with AD and partially with ALS.

## Two distinct pathways connect AD genetics to peripheral organs

Integrating the cell-type and gene-level analyses, we identified two mechanistically distinct pathways through which AD genetic risk operates in peripheral organs (**Figure 6**):

1. **Myeloid/immune pathway** (TREM2, CD33, SPI1, INPP5D, PLCG2): Active in all organs through tissue-resident macrophages. Brain PCC = 0.93–0.97; Lung = 0.59–0.87; Liver = 0.59–0.88. Functionally linked to phagocytosis, innate immunity, and Aβ clearance.

2. **APP processing pathway** (APP, SORL1, PSEN1, ADAM10, CD2AP, CR1, PICALM): Selectively active in colonic epithelium. Colon PCC = 0.70–0.82; Brain = 0.19–0.94 (variable by gene). Functionally linked to amyloid precursor protein processing, endosomal trafficking, and potentially peripheral Aβ production.

These two gene sets showed distinct organ-specificity profiles: myeloid genes were highest in Brain/Lung/Liver/Heart, while APP pathway genes were highest in Colon. The convergence of both pathways in brain (where microglia express both gene sets) may explain why AD manifests primarily as a brain disease despite having systemic genetic underpinnings.
