# Paper Outline

## Story Arc
**Background**: AD GWAS identified 75+ risk loci enriched in immune/myeloid pathways, but prior studies used bulk/single-cell data without spatial resolution or whole-body perspective.
**Gap**: No study has systematically mapped AD genetic risk across all organs at spatial resolution, nor examined how this pattern changes with aging or differs from other neurodegenerative diseases.
**Approach**: gsMap + whole-body Array-seq + multi-age CMAP mapping + multi-trait comparison.
**Key Findings**:
1. AD risk is myeloid-specific across ALL 15 organs (not just brain)
2. Gut epithelial cells show unexpected AD enrichment independent of immune cells
3. This gut pattern is AD-specific — absent in PD, partial in ALS
4. Pattern stable from 1-30 months (no age dependence)
**Implication**: Host genetic evidence for the gut-brain axis in AD; AD has unique peripheral risk mechanisms.

---

## Figure Plan (7 main figures for Cell)

### Figure 1: Whole-Body AD Risk Atlas
- **1A**: Study design schematic (Array-seq → gsMap → CMAP → multi-trait)
- **1B**: Whole-body spatial map showing all 15 organs (from `Wholebody_composite_young.png` — right panel)
- **1C**: Organ-level enrichment bar chart (from `Fig3A_organ_barplot_cauchy.png`)
- **Message**: AD genetic risk is distributed across the entire body, strongest in immune organs

### Figure 2: Myeloid Cells as Universal Carriers of AD Risk
- **2A**: Cell type × organ dot plot (from `Fig4A_dotplot_celltype_organ.png`)
- **2B**: Per-organ: immune vs non-immune enrichment comparison (new figure needed — grouped bar chart)
- **2C**: Tissue-resident macrophage comparison: Kupffer (Liver) vs Alveolar (Lung) vs Microglia (Brain) — all highly enriched
- **Message**: In every organ, myeloid/macrophage cells carry the strongest AD risk — this is a lineage-specific, not organ-specific, phenomenon

### Figure 3: Brain Microglia Enrichment
- **3A**: Brain cell type bar chart — microglia p<10⁻⁸, all neurons p>0.85 (from `AD_vs_PD_brain_celltypes.png` left panel)
- **3B**: Brain spatial risk map across 6 ages with MG p labels (from `Age_progression_Brain.png`)
- **3C**: Brain vs immune organs trajectory (from `Brain_microglia_vs_nonmyeloid.png` left panel)
- **Message**: Brain's AD enrichment is entirely microglia-driven; after enrichment, brain matches peripheral immune organs

### Figure 4: Gut Epithelial Enrichment — A Novel AD-Specific Finding
- **4A**: Colon/Ileum/Stomach: immune vs epithelial enrichment (AD vs PD) (from `AD_vs_PD_gut_celltype.png`)
- **4B**: Top cell types in gut organs — AD shows secretory cell, enterocyte enrichment
- **4C**: Known AD risk genes expressed in gut epithelial cells (new analysis: which specific genes?)
- **Message**: Gut epithelial cells independently enrich for AD risk variants — a non-immune pathway supporting the gut-brain axis hypothesis

### Figure 5: Age-Resolved Stability
- **5A**: Organ × age heatmap (from `Age_AD_risk_heatmap_final.png`)
- **5B**: Whole-body spatial maps at 3m vs 21m (from `Wholebody_composite_young_vs_old.png`)
- **5C**: Spearman correlation statistics — all organs flat
- **Message**: AD genetic risk landscape is age-independent, consistent with germline-encoded susceptibility

### Figure 6: Multi-Trait Comparison — AD vs PD vs ALS
- **6A**: Three-disease organ-level comparison (from `AD_PD_ALS_organ_comparison.png`)
- **6B**: Myeloid enrichment gradient heatmap: AD >> ALS > PD (from `Myeloid_enrichment_AD_ALS_PD.png`)
- **6C**: Brain cell types: AD (microglia) vs PD (oligodendrocyte) (from `AD_vs_PD_brain_celltypes.png`)
- **Message**: Myeloid enrichment strength and gut epithelial involvement are disease-specific features, not general neurodegeneration

### Figure 7: Integrated Model
- **7A**: Schematic: DNA variants → myeloid dysfunction → 3 brain pathways + gut epithelial pathway
- **7B**: Summary heatmap: organ × trait × cell type
- **Message**: AD genetic architecture operates through two distinct peripheral pathways — myeloid and gut epithelial — with disease-specific patterns

---

## Table Plan

### Table 1: Data Sources and Coverage
- 8 scRNA-seq sources, organs, ages, cell counts, GEO accessions

### Table 2: Multi-Trait Organ Enrichment Summary
- 15 organs × 5 traits (AD/PD/ALS/MS/FTD) — min Cauchy p

---

## STAR Methods Subsections
1. Spatial transcriptomics data (Array-seq, Cell 2026)
2. GWAS summary statistics (AD, PD, ALS, MS, FTD)
3. Single-cell RNA-seq data sources (PanSci, TMS, MCA2, GSE179358, etc.)
4. gsMap pipeline (quick_mode, parameters)
5. CMAP spatial mapping (run_cmap_fast, parameters)
6. Brain microglia enrichment strategy
7. Age-resolved analysis pipeline
8. Multi-trait comparison pipeline
9. Statistical analysis (Cauchy combination, Spearman correlation)
10. Software and data availability

---

## Supplementary Plan
- **Fig S1**: QC violins, MT%, replicate consistency
- **Fig S2**: gsMap validation (GVAE UMAP, QQ plot)
- **Fig S3**: Per-organ spatial risk maps (16 organs)
- **Fig S4**: All 15 organ age progression panels
- **Fig S5**: PD and ALS organ × age heatmaps
- **Fig S6**: MS and FTD results (when available)
- **Fig S7**: Cross-species homolog coverage
- **Table S1**: Full 113 organ × age results for AD
- **Table S2**: Cell type enrichment per organ (all annotations)
- **Table S3**: Brain cell type decomposition per age
