# Notebooks — analysis pipeline + figures

The notebooks below preserve the **full analysis pipeline** in execution order. Notebooks tagged **Figure** ship with pre-rendered output cells so you can read them top-to-bottom without execution. Notebooks tagged **Pipeline** drive long-running computations (gsMap S-LDSC, CMAP projection, full-summary-statistics coloc) and need to be re-executed against the data sources listed in `../data_manifest/README.md`.

Removing any of these notebooks would break the analysis logic — keep the full set.

| Notebook | role | code | figures | what it does |
|---|---|---:|---:|---|
| **Phase 1 — Data** | | | | |
| `01_data_preprocessing.ipynb` | Pipeline | 5 | *pipeline-only* | 01 Data Preprocessing |
| `02_gsmap_analysis.ipynb` | Pipeline | 3 | *pipeline-only* | 02 gsMap Analysis: AD GWAS Projection onto Whole-Body Spatial Transcriptomics |
| `03_cmap_age_mapping.ipynb` | Pipeline | 4 | *pipeline-only* | 03 CMAP Age-Resolved Spatial Mapping |
| **Phase 2 — Figures** | | | | |
| `Paper1_Fig1_wholebody_atlas.ipynb` | Figure | 10 | **8** | Paper 1 — Figure 1: Whole-Body AD Genetic Risk Atlas |
| `Paper1_Fig2_colon_exception.ipynb` | Figure | 10 | **8** | Paper 1 — Figure 2: Colon Epithelium is the Sole Exception to Myeloid Dominance |
| `Paper1_Fig4_cross_validation.ipynb` | Figure | 8 | **7** | Figure 4: Cross-Species Cross-Dataset Validation + Barrier Mechanism |
| `Paper1_Fig5_multi_trait.ipynb` | Figure | 9 | **8** | Figure 5: 8-Trait Multi-Disease ComparisonCompare gsMap spatial enrichment acros |
| **Phase 3 — Genetics** | | | | |
| `21_coloc_full_eqtl.ipynb` | Pipeline | 10 | *pipeline-only* | 21. Colocalization Analysis: AD GWAS × Colon eQTL (Full Summary Statistics)## Pu |
| `21_eqtl_coloc_analysis.ipynb` | Figure | 8 | **2** | 21. eQTL & Colocalization Analysis: AD GWAS x Human Colon |
| `21_eqtl_coloc_analysis_executed.ipynb` | Figure | 8 | **2** | 21. eQTL & Colocalization Analysis: AD GWAS x Human Colon |
| **Phase 4 — Supp** | | | | |
| `S01_qc.ipynb` | Figure | 7 | **4** | Supplementary Figure S1: Data Quality Control |
| `S02_gsmap_validation.ipynb` | Figure | 7 | **4** | Supplementary Figure S2: gsMap Technical Validation |
| `S03_per_organ_maps.ipynb` | Figure | 4 | **1** | Supplementary Figure S3: Per-Organ Detailed Spatial Risk Maps |
| `S04_extended_celltype.ipynb` | Figure | 6 | **3** | Supplementary Figure S4: Extended Cell-Type Analysis |
| `S05_extended_genes.ipynb` | Figure | 6 | **2** | Supplementary Figure S5: AD Risk Gene Extended Analysis |
| `S06_sensitivity.ipynb` | Figure | 6 | **3** | Supplementary Figure S6: Sensitivity Analysis |
| `S07_cross_species.ipynb` | Figure | 7 | **3** | Supplementary Figure S7: Cross-Species Validation |
| **Phase 5 — Master** | | | | |
| `ng_fig_gsmap_eqtl_verify.ipynb` | Figure | 36 | **19** | Whole-Body Spatial Mapping of AD Genetic Risk — Nature Genetics Figures (v2) |
| `paper_figure.ipynb` | Figure | 43 | **25** | Paper figure |
| `paper_figure_supp.ipynb` | Figure | 2 | **1** | Supplementary / extended-data figures (Para 3 S3: organ-wise Young vs Old effect size with CI + sample-level variability) |

## Re-running pipeline notebooks

Pipeline notebooks (01, 02, 03, 21_coloc_full_eqtl) need data not bundled here. Get the data per `../data_manifest/README.md`, then:

```bash
pip install -r ../requirements.txt
jupyter nbconvert --to notebook --execute --inplace 01_data_preprocessing.ipynb
# (or just open in Jupyter and run all cells)
```

## Re-rendering figure notebooks

Figure notebooks read the **intermediate outputs** produced by the pipeline notebooks (CSVs in `../results/`, h5ad in `../data/age_merged/` once you build it). Re-render with the same `nbconvert --execute` command above.
