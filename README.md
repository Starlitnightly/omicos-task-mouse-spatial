# omicos-task-mouse-spatial

Reproducible analysis code for the paper
**"Whole-body spatial mapping of Alzheimer's disease genetic risk reveals a colonic epithelial APP-processing pathway"**

> One-line abstract — Spatially-resolved partitioning of AD GWAS heritability across 16 mouse organs uncovers two genetically distinct programs: a universal tissue-resident myeloid axis (TREM2/CD33/SPI1/INPP5D/PLCG2) and a colon-exclusive epithelial APP-processing axis (APP/SORL1/PSEN1/ADAM10/CD2AP/CR1/PICALM), with GTEx colon eQTL coloc + SharePro + SMR nominating CR1, PICALM, CD2AP, and ADAM10 as colon-intrinsic causal candidates.

---

## What's in this repo

| Path | Purpose |
|------|---------|
| `notebooks/` | End-to-end Jupyter notebooks for preprocessing, gsMap, eQTL colocalization, multi-trait comparison, and every main/supplementary figure. |
| `scripts/` | Headless Python pipelines used on the cluster — aggregation, h5ad merging, per-organ processing, and figure-generation scripts that mirror the notebooks. |
| `figures/` | Output directory for paper figures (PNG/PDF/SVG) produced by the notebooks and scripts. |
| `results/` | All derived tables: Cauchy summaries per organ × trait, spot-level p-values, coloc / SharePro / SMR outputs, GO enrichments. |
| `docs/` | Scientific narrative (`AD_colon_story.md`), conversation walkthrough, and paper drafts. |
| `data_manifest/` | Pointers to public datasets used (raw data is **not** redistributed). |

---

## Two main results

### 1. Universal myeloid axis (the rule)

Across 15 mouse organs, tissue-resident myeloid populations carry AD heritability:

- **11 of 15 organs** have a macrophage/monocyte annotation as the top cell type (Liver Kupffer p = 4.9 × 10⁻¹⁰, Lung alveolar Mφ p = 1.2 × 10⁻¹⁰, BM monocyte p = 5.2 × 10⁻¹⁰, brain microglia p = 2.4 × 10⁻⁸).
- All neuronal subtypes are non-significant (p > 0.85).
- Driver genes (per-spot PCC): **TREM2, CD33, SPI1, INPP5D, PLCG2**, peaking in brain (PCC 0.93–0.97), lung, and liver.

### 2. Colon-only epithelial APP-processing pathway (the exception)

Colon (and to a lesser extent ileum/stomach) flips this pattern entirely:

- Top cell types are **colon secretory** (p = 7.3 × 10⁻⁷) and **enterocytes** (p = 3.8 × 10⁻⁶) — non-immune.
- Driver genes form a distinct **APP-processing module**: **APP, SORL1, PSEN1, ADAM10, CD2AP, CR1, PICALM** (PCC 0.70–0.82 in colon).
- Colon is the **only** organ where the APP module's mean PCC (0.68) exceeds the myeloid module's (0.27).
- GTEx colon eQTL fine-mapping → **CR1** (coloc PP.H4 = 0.993, sigmoid), **PICALM** (SharePro share = 0.998), **CD2AP** (SharePro), **ADAM10** (SMR p = 4 × 10⁻⁴) share causal variants with AD GWAS in colon tissue.
- Replicated across Burclaff/Smillie scRNA, Oliveira human colon Visium, Das mouse colon Visium, Stereo-seq aging atlas, and Xenium colon panels.
- Specificity controls: cell-composition downsampling does not erase the signal; multi-trait scan shows AD ≫ ALS > PD > others in colonic epithelium; effect is age-stable (1–30 months).

Full evidence chain: [`docs/AD_colon_story.md`](docs/AD_colon_story.md).

---

## Repository structure

```
omicos-task-mouse-spatial/
├── README.md
├── notebooks/
│   ├── 01_data_preprocessing.ipynb         # Build h5ad inputs for gsMap
│   ├── 02_gsmap_analysis.ipynb             # Run gsmap quick_mode per organ × replicate
│   ├── 03_cmap_age_mapping.ipynb           # Age-resolved (1–30 mo) gsMap projection
│   ├── 21_coloc_full_eqtl.ipynb            # GTEx colon eQTL × Bellenguez 2022 coloc.abf
│   ├── 21_eqtl_coloc_analysis.ipynb        # Coloc + SharePro + SMR pipeline
│   ├── ng_fig_gsmap_eqtl_verify.ipynb      # Cross-dataset validation (Visium, sc, Stereo-seq)
│   ├── Paper1_Fig1_wholebody_atlas.ipynb   # Fig 1 — whole-body heritability map
│   ├── Paper1_Fig2_colon_exception.ipynb   # Fig 2 — colon flips the myeloid rule
│   ├── Paper1_Fig4_cross_validation.ipynb  # Fig 4/5 — coloc/SharePro/SMR + sc atlas
│   ├── Paper1_Fig5_multi_trait.ipynb       # Multi-trait AD vs ALS/PD/T2D/SCZ/EA/height
│   ├── paper_figure.ipynb                  # Consolidated figure builder
│   └── S0[1-7]_*.ipynb                     # Supplementary: QC, validation, per-organ maps,
│                                           # extended celltype/gene panels, sensitivity, cross-species
├── scripts/
│   ├── aggregate_results.py                # Roll up per-organ gsMap → Cauchy tables
│   ├── aggregate_age_results.py            # Same for age-resolved CMAP samples
│   ├── build_age_provenance*.py            # Sample-provenance ledgers for age cohort
│   ├── merge_*_h5ad.py                     # Stitch raw spatial h5ads into per-organ inputs
│   ├── preprocess_new_sc.py                # Human colon scRNA preprocessing (Burclaff/Smillie)
│   ├── process_brain_microglia_enriched.py # Microglia-enriched re-annotation for brain
│   ├── process_mca2.py                     # Mouse Cell Atlas v2 annotation harmonization
│   ├── process_remaining_gaps.py           # Gap-fill for organs missing initial annotations
│   ├── figure1_paper.py … figure7.py       # Headless figure builders (mirror notebooks)
│   └── figS[1-8]*.py                       # Supplementary figure builders
├── figures/                                # Output figures (populated by notebooks/scripts)
├── results/
│   ├── all_organs_cauchy_AD.csv.gz         # Per organ × annotation Cauchy combined p
│   ├── all_organs_cauchy_{ALS,PD,MS,FTD}.csv.gz   # Multi-trait Cauchy tables
│   ├── all_organs_spot_AD_pvalues.csv.gz   # Per-spot −log10(p) for AD
│   ├── age_all_cauchy_AD.csv.gz            # Age-resolved Cauchy (126 samples)
│   ├── age_all_spots_AD_pvalues.csv.gz     # Age-resolved spot p-values
│   ├── age_sample_provenance*.csv          # Sample → source ledger for age cohort
│   ├── all_8traits_organ_best.csv          # Best-organ summary across 8 GWAS traits
│   ├── coloc_FULL_eqtl_results.csv         # coloc.abf for 10 AD locus genes × GTEx colon
│   ├── sharepro/                           # Per-gene SharePro credible-set outputs
│   ├── smr_colon_AD_results.csv            # SMR + HEIDI for colon eQTL × AD
│   ├── gtex_eqtl_query_results.csv         # Raw GTEx eQTL queries
│   └── GO_enrichment_colon_AD_top200.csv   # GO enrichment of top colon AD-PCC genes
├── docs/
│   ├── AD_colon_story.md                   # Full scientific narrative (figure ↔ evidence map)
│   ├── conversation_walkthrough.md         # Step-by-step analytical reasoning
│   └── paper/                              # Manuscript sections (IMRAD)
└── data_manifest/                          # Public data accession pointers
```

---

## How to run

### Environment

Core dependencies (Python 3.10+):

```bash
conda create -n omicos-mouse-spatial python=3.10 -y
conda activate omicos-mouse-spatial

# spatial + single cell
pip install scanpy anndata squidpy

# gsMap (https://github.com/JianYang-Lab/gsMap)
pip install gsMap

# colocalization / fine-mapping
# - coloc (R package) via Bioconductor or CRAN
# - SharePro: https://github.com/zhwm/SharePro_coloc
# - SMR: https://yanglab.westlake.edu.cn/software/smr/

# stats / plotting
pip install numpy pandas scipy statsmodels matplotlib seaborn scikit-learn
```

R packages: `coloc`, `data.table`, `dplyr`, `ggplot2`.

External binaries: `smr_Linux` on `$PATH`; `magma` v1.10; `plink` 1.9 / 2.

### Step-by-step

1. **Fetch data.** Follow `data_manifest/` to download the public datasets listed below into a local `data/` tree (not committed). Expected layout matches the paths in `docs/AD_colon_story.md` § 8.
2. **Preprocess spatial inputs.** `notebooks/01_data_preprocessing.ipynb` (or `scripts/merge_original_h5ad.py` + `merge_age_h5ad.py`).
3. **Run gsMap per organ × replicate.** `notebooks/02_gsmap_analysis.ipynb` → outputs under `models/gsmap_output/<organ>/<replicate>/`. Age-resolved cohort uses `notebooks/03_cmap_age_mapping.ipynb`.
4. **Aggregate to Cauchy / spot tables.** `python scripts/aggregate_results.py` and `python scripts/aggregate_age_results.py` → populates `results/all_organs_cauchy_*.csv.gz` and `results/age_all_*.csv.gz`.
5. **Colocalization on AD locus genes.** `notebooks/21_coloc_full_eqtl.ipynb` (coloc.abf) → `results/coloc_FULL_eqtl_results.csv`. Then SharePro and SMR — see `notebooks/21_eqtl_coloc_analysis.ipynb` → `results/sharepro/` + `results/smr_colon_AD_results.csv`.
6. **Cross-dataset validation.** `notebooks/ng_fig_gsmap_eqtl_verify.ipynb` runs through Burclaff/Smillie sc, Oliveira/Das Visium, Stereo-seq aging, and Xenium panels.
7. **Build figures.** Run `notebooks/Paper1_Fig*.ipynb` and `notebooks/S0*.ipynb`, or use the headless `scripts/figure*.py` / `scripts/figS*.py` equivalents. Outputs land in `figures/`.

A full run (steps 3–7) takes ~24–48 h on a single GPU node; gsMap step dominates wall-clock.

---

## Data sources

Raw data is **not** included in this repository. All datasets are publicly available — see `data_manifest/` for exact accession URLs and per-file SHAs.

| Dataset | Accession / Source | Used for |
|---|---|---|
| Whole-body mouse Array-seq (Clevenger et al. 2026) | **GEO GSE248904** | Primary spatial input — 2 sagittal C57BL/6J sections, ~1.2M spots, 16 organs |
| AD GWAS (Bellenguez et al. 2022, *Nat Genet*) | EBI GWAS Catalog **GCST90027158** | AD heritability partitioning (cases ≈ 90k, total ≈ 788k) |
| Additional GWAS (PD, ALS, MS, FTD, T2D, SCZ, height, EA) | GWAS Catalog / IGAP / PGC / GIANT / SSGAC | Multi-trait specificity comparison |
| GTEx v8 (Colon — Transverse + Sigmoid) | **GTEx Portal** (eQTL all-pairs, fine-mapping) | coloc / SharePro / SMR causal-variant analysis |
| Human colon single-cell (Burclaff et al. 2022) | **GSE185224** | Epithelial-enriched scRNA validation |
| Human colon single-cell (Smillie et al. 2019) | **SCP259 / GSE121380** | Immune compartment scRNA control |
| Human colon Visium (Oliveira et al. 2024, P3NAT) | **Zenodo / GSE-listed** | Human spot-level epithelial coupling (Fig 3E) |
| Mouse colon Visium (Das et al. 2022) | **GSE190037** | Cross-species spatial validation |
| Stereo-seq mouse aging atlas | **CNGB CNP0002671** / **GSE247719** (proxy) | Age-stability of APP module in intestine vs hippocampus |
| Human colon Xenium | 10x Genomics public Xenium colon panel | Morphological anchor (APOE-only relevant probe; CR1/PICALM/CD2AP/ADAM10 not on panel — documented caveat) |
| Tabula Muris Senis | **FigShare 8273102** | Mouse single-cell reference for annotation harmonization |
| Mouse Cell Atlas v2 | **figshare / ABA** | Per-organ celltype annotation backbone |

For convenience, `data_manifest/` includes a small JSON registry of each dataset's accession, expected local path, and the notebooks that consume it.

---

## Citations

If you use this code or its derived tables, please cite the paper above and the underlying data/tools:

- gsMap — Song et al. 2024 (*Nat Genet*)
- S-LDSC / MAGMA — Finucane et al. 2015; de Leeuw et al. 2015
- Bellenguez C, et al. New insights into the genetic etiology of Alzheimer's disease. *Nat Genet* 54, 412–436 (2022).
- coloc — Giambartolomei et al. 2014; Wallace 2020
- SharePro — Zhang & Greenwood 2023
- SMR / HEIDI — Zhu et al. 2016 (*Nat Genet*)
- GTEx Consortium 2020 (*Science*)
- Clevenger et al. 2026 — whole-body Array-seq (GSE248904)
- Burclaff et al. 2022; Smillie et al. 2019; Oliveira et al. 2024; Das et al. 2022; Stereo-seq aging atlas (CNGB / GSE247719)
- Tabula Muris Senis Consortium 2020

A consolidated BibTeX file lives at `docs/paper/references.bib`.

---

## License

[MIT](LICENSE) — code only. Public datasets retain their original licenses; consult each repository's terms of use.

---

## Acknowledgments

We used **omicos** (omicos, Opus 4.x) as an analytical assistant for prototyping notebooks, debugging pipelines, and drafting figure code. All scientific claims and final code were reviewed and validated by the authors. Compute provided by Stanford Sherlock.
