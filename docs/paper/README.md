# Whole-Body Spatial Mapping of AD Genetic Risk

## Paper Info

| Item | Value |
|------|-------|
| **Working Title** | Whole-body spatial mapping reveals Alzheimer's disease genetic risk concentrated in myeloid cells with disease-specific gut epithelial enrichment |
| **Paper Type** | Original Article |
| **Target Journal** | Cell |
| **Language** | English |
| **Research Question** | Using gsMap to project AD GWAS onto whole-body spatial transcriptomics, we discovered that AD genetic risk is concentrated in myeloid lineage cells across all organs, with a novel gut epithelial enrichment not seen in PD — providing spatial genetic evidence for the gut-brain axis in AD. |
| **Reporting Guideline** | STROBE (observational/computational genomics) |
| **Created** | 2026-04-10 |
| **Last Updated** | 2026-04-10 |

## Journal Requirements

| Item | Value |
|------|-------|
| Word Limit (total) | 7,000 (main text + figure legends, excluding STAR Methods) |
| Abstract Word Limit | 150 words (Summary, unstructured) |
| Abstract Format | Unstructured ("Summary" in Cell style) |
| Citation Style | Cell Press (numbered, superscript in text) |
| Figure/Table Limit | 7 main figures/tables |
| Required Sections | Summary, Introduction, Results, Discussion, STAR Methods, Supplemental Information |
| Graphical Abstract | Required |
| Keywords | Not required (Cell uses "Highlights" instead — 3-4 bullet points, max 85 characters each) |
| AI Disclosure Required | Yes, in STAR Methods or Acknowledgments |
| Special Requirements | eTOC blurb (max 80 words), Highlights (3-4 bullets), Lead Contact, Data/Code Availability |

## Project Status

| Phase | Status | Last Updated | Notes |
|-------|--------|-------------|-------|
| Ethics & Protocol | N/A | - | Computational study, public data only |
| Data Organization | Done | 2026-04-10 | 8 data sources, all public GEO/figshare |
| Literature Search | Not Started | - | `00_literature/` |
| Outline | Not Started | - | `01_outline.md` |
| Tables & Figures | In Progress | 2026-04-10 | All analysis figures generated in `../figures/` |
| Methods & Results | Not Started | - | `sections/02_methods.md`, `sections/03_results.md` |
| Introduction & Conclusion | Not Started | - | `sections/04_introduction.md`, `sections/06_conclusion.md` |
| Discussion | Not Started | - | `sections/05_discussion.md` |
| Abstract & Title | Not Started | - | `sections/07_abstract.md`, `sections/08_title.md` |
| Humanize | Not Started | - | Phase 4 |
| References | Not Started | - | `references/09_references.md` |
| Quality Review | Not Started | - | `checklists/` |
| Co-author Review | Not Started | - | `coauthor-review/` |
| Pre-Submission | Not Started | - | `submissions/v1_cell/` |

## Data Management

| Item | Value |
|------|-------|
| Raw data location | `./data/  (not in repo; see data_manifest/)` |
| Data format | h5ad (AnnData), sumstats.gz (GWAS), csv.gz (results) |
| Analysis software | Python 3.10, gsMap 1.73.7, CMAP_py 1.1.0, omicverse 2.1.1 |
| De-identified | N/A (mouse data, public) |

### Key Figures (generated)
- `Age_AD_risk_heatmap_final.png` — organ × age enrichment heatmap
- `Age_AD_risk_trajectory_final.png` — age trajectory all organs
- `Brain_microglia_vs_nonmyeloid.png` — brain cell type comparison
- `Age_progression_Brain.png` — brain across ages
- `Wholebody_composite_young.png` — whole-body ST + CMAP + risk
- `Wholebody_per_age_composite_final.png` — per-age whole-body
- `AD_PD_ALS_organ_comparison.png` — multi-trait organ comparison
- `Myeloid_enrichment_AD_ALS_PD.png` — myeloid gradient
- `AD_vs_PD_gut_celltype.png` — gut immune vs epithelial
- `AD_vs_PD_brain_celltypes.png` — brain cell types AD vs PD
- 15× `Age_progression_{organ}.png` — per-organ progressions

### Key Results
- `results/age_all_cauchy_AD.csv.gz` — 126 AD Cauchy results
- `results/age_organ_summary.csv` — organ × age summary

## Key Decisions Log

| Date | Decision | Rationale |
|------|----------|-----------|
| 2026-04-10 | Project initialized | All analysis complete, ready for manuscript |
| 2026-04-10 | Target Cell | Whole-body spatial + multi-trait + novel gut finding = high impact |
| 2026-04-10 | Focus on gut epithelial as novel finding | Data-driven: Colon/Ileum epithelial enrichment is AD-specific vs PD/ALS |
