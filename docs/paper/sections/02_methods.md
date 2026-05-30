# STAR Methods

## Resource Availability

### Lead Contact
Further information and requests for resources should be directed to the lead contact.

### Data and Code Availability
- Whole-body spatial transcriptomics data: GEO GSE248904
- AD GWAS summary statistics: Bellenguez et al. 2022 (GCST90027158)
- PD GWAS: Nalls et al. 2019 (GCST009325)
- ALS GWAS: van Rheenen et al. 2021 (GCST90027163)
- Single-cell RNA-seq data: PanSci (GSE247719), TMS (GSE132042), MCA2 (GSE153562), hippocampal microglia (GSE179358), aging colon (GSE168448), Calico spleen (GSE132901)
- gsMap software: https://github.com/LeonSong1995/gsMap (v1.73.7)
- Analysis code: [to be deposited]

## Method Details

### Whole-body spatial transcriptomics
We used the Array-seq whole-body mouse spatial transcriptomics atlas (GSE248904), which profiled sagittal sections of 6-week-old C57BL/6J mice using DNA microarray-based spatial capture. The dataset contains two biological replicates (CTRL1 and CTRL2), each covering 16 organ systems with ~600,000 spots per section. Each spot was annotated with organ identity and sub-regional annotation. Gene expression was stored as raw UMI counts in the `count` layer. We used the `x_plotting` and `y_plotting` coordinates for whole-body spatial visualization, and per-organ h5ad files with local `spatial` coordinates for gsMap analysis.

### GWAS summary statistics
We obtained GWAS summary statistics for five traits:
- **Alzheimer's disease (AD)**: Bellenguez et al. 2022, 111,326 cases and 677,663 controls, 75 genome-wide significant loci
- **Parkinson's disease (PD)**: Nalls et al. 2019, up to 482,730 subjects
- **ALS**: van Rheenen et al. 2021, 29,612 cases and 122,656 controls
- **Multiple sclerosis (MS)**: FinnGen R12, 2,926 cases and 495,931 controls
- **Frontotemporal dementia (FTD)**: Ferrari et al. 2014, 2,154 cases and 4,308 controls

All summary statistics were formatted to the gsMap-required format (SNP, A1, A2, Z, N) using standard conversion: Z = beta/SE for studies reporting beta and standard error, or Z derived from log(OR)/SE for studies reporting odds ratios.

### gsMap spatial GWAS projection
We applied gsMap (v1.73.7) in quick_mode to project GWAS signals onto spatial transcriptomics data. gsMap assesses whether SNPs located near genes highly expressed at a given spatial spot are enriched for trait heritability, using the stratified LD score regression (S-LDSC) framework. For each spot, gsMap: (1) learns latent representations using a graph attention autoencoder that captures spatial gene expression patterns; (2) computes Gene Specificity Scores (GSS) reflecting each gene's spatial expression specificity; (3) generates per-spot LD score annotations using pre-computed SNP-gene weight matrices; (4) performs S-LDSC regression to test heritability enrichment; and (5) applies Cauchy combination to aggregate spot-level p-values within each annotated cell type or region. Mouse-to-human gene mapping used the homolog file provided with gsMap resources.

For baseline analysis, we ran gsMap on all 32 per-organ h5ad files (16 organs × 2 replicates) for each of the five GWAS traits. For age-resolved analysis, gsMap was run on CMAP-mapped per-age h5ad files (see below).

### Multi-age single-cell data integration
To extend the analysis across mouse aging, we integrated single-cell RNA-seq data from eight sources covering 15 organs and 10 age points (1–30 months):

| Source | Organs | Ages | Method | GEO |
|--------|--------|------|--------|-----|
| PanSci | 9 non-immune organs + Brain | 3/6/12/16/23m | EasySci snRNA-seq | GSE247719, GSE212606 |
| TMS FACS | 23 organs | 3/18/24m | FACS scRNA-seq | GSE132042 |
| TMS Droplet | 16 organs | 1/3/18/21/24/30m | 10X Chromium | GSE132042 |
| MCA2 | Stomach/Intestine/Spleen/Thymus | 12/18/24m | Microwell-seq | GSE153562 |
| GSE179358 | Brain (hippocampal microglia) | 6/12/18/24m | 10X Chromium | GSE179358 |
| GSE168448 | Colon (aging epithelial) | ~21m | 10X Chromium | GSE168448 |
| GSE132901 | Spleen | 7m (~6m proxy) | 10X Chromium | GSE132901 |

For each organ and age, we subsampled up to 10,000 cells, filtered to WT/control conditions, and subset to genes shared with the corresponding spatial transcriptomics data.

### CMAP spatial mapping
We used CMAP_py (v1.1.0) to map age-resolved single cells onto the 6-week Array-seq spatial coordinates. CMAP performs hierarchical spatial mapping: (1) Harmony integration of single-cell and spatial data in PCA space; (2) SVM-based domain classification to assign cells to spatial annotations; (3) within-domain nearest-neighbor positioning to assign precise spatial coordinates. Parameters: `spatial_data_type='square'`, `svm_threshold=0.5`, `num_epochs=2000`, `n_near_spot=5`, `radius=1/6`. For organs with >25,000 spots, spatial data was subsampled to 15,000–25,000 spots. Mapped cells with valid spatial coordinates were retained and saved as per-age h5ad files with `obsm['spatial']` and `layers['count']`.

### Brain microglia enrichment
PanSci Brain data (EasySci snRNA-seq) contained only ~1.2% microglia among total cells, resulting in weak AD enrichment due to signal dilution. To address this, we enriched each Brain age point with FACS-sorted microglia from TMS Brain_Myeloid (3/18/24m) and hippocampal CD11b+ cells from GSE179358 (6/12/18/24m). For each age, we combined: (1) PanSci non-microglia cells (~5,000, providing spatial neuronal/glial context); (2) TMS FACS microglia (~1,500–3,000, providing strong myeloid signal); (3) GSE179358 hippocampal microglia (~1,500–3,000). For the 21-month time point, where neither TMS nor GSE179358 had data, we used TMS 18m and 24m microglia as a proxy combined with PanSci 21m full microglia complement (3,152 cells from the complete dataset). This enrichment strategy achieved 30–60% microglia content per age point while maintaining diverse brain cell type representation.

### Gene-level analysis
gsMap produces per-gene diagnostic information including Pearson Correlation Coefficient (PCC) between each gene's spatial marker score and the GWAS trait association signal. We used these PCC values to identify which AD risk genes drive enrichment in each organ. Known AD GWAS genes were classified into two categories based on their cross-organ PCC patterns: (1) myeloid pathway genes (TREM2, CD33, SPI1, INPP5D, ABI3, PLCG2, MS4A6A, MEF2C, BIN1) — those with highest PCC in Brain, Lung, Liver, and Heart; and (2) APP processing pathway genes (APP, SORL1, PSEN1, ADAM10, CD2AP, PICALM, CR1) — those with uniquely high PCC in Colon.

### Statistical analysis
For each organ and annotation, significance was assessed using the Cauchy combination test, which aggregates per-spot p-values while accounting for correlation structure. Age trends were evaluated using Spearman rank correlation of -log10(Cauchy p) versus age (months) for each organ, with p < 0.05 considered significant. Multi-trait comparisons used the minimum Cauchy p-value across all ages for each organ-trait combination. All analyses were performed using Python 3.10 with gsMap 1.73.7, scanpy 1.10, and anndata 0.10.
