# Data Manifest — Whole-Body Mouse Spatial × AD GWAS Project

This manifest lists all external datasets required to reproduce the analyses in
the `omicos-task-mouse-spatial` project. Every entry contains the canonical
source citation, accession, file count, approximate total size, and a direct
download URL. Where applicable, the resolved snapshot date is noted; for live
APIs (e.g. GTEx) the query parameters are given instead of a file dump.

Total approximate download volume: ~75 GB raw + ~30 GB derived bundles.

---

## 1. Primary spatial transcriptomics

### 1.1 Whole-body mouse Array-seq (anchor dataset)
- **Name:** Whole-body Array-seq mouse atlas
- **Paper:** Clevenger J. et al., *Cell* 2026 (in press) — "Whole-body spatial transcriptomics of the adult mouse via Array-seq"
- **Accession:** GEO `GSE248904`
- **Files:** 126 spatial sections (Visium-format `filtered_feature_bc_matrix.h5` + `tissue_positions_list.csv` + image), ~252 files total
- **Size:** ~3.4 GB compressed
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248904
- **Notes:** This is the anchor dataset used as the gsMap spatial input for all 126 sections.

### 1.2 Stereo-seq mouse aging atlas
- **Name:** Mouse Aging Stereo-seq atlas (MOSTA-aging)
- **Paper:** Chen A. et al., *Cell* 2024 — "Spatiotemporal transcriptomic atlas of mouse aging"
- **Accession:** CNGB `STDS0000247`
- **Files:** ~40 Stereo-seq slides (`.h5ad` + `.gef`)
- **Size:** ~18 GB
- **URL:** https://db.cngb.org/stomics/project/STDS0000247
- **Notes:** Used for cross-platform validation of microglial/gut signatures across age.

### 1.3 Xenium 5K human colon
- **Name:** Xenium Human Colon Preview (5K panel)
- **Paper:** 10x Genomics application note, 2024
- **Accession:** 10x Genomics public dataset (no GEO)
- **Files:** 1 Xenium output bundle (`cells.parquet`, `transcripts.parquet`, `cells.zarr.zip`, morphology image)
- **Size:** ~6 GB
- **URL:** https://www.10xgenomics.com/datasets/xenium-prime-5k-human-colon-ffpe
- **Notes:** Human cross-species validation of the gut APP-processing pathway.

### 1.4 Das 2022 mouse Visium colon
- **Name:** Mouse colon Visium (inflammation series)
- **Paper:** Das S. et al., *Nat Commun* 2022 — "Spatial transcriptomics of healthy and inflamed mouse colon"
- **Accession:** GEO `GSE189184`
- **Files:** 8 Visium sections
- **Size:** ~2.1 GB
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189184

### 1.5 Oliveira P3NAT human Visium
- **Name:** P3NAT human gut Visium atlas
- **Paper:** Oliveira M.F. et al., *Nature* 2024 — "Characterization of immune cell populations in the tumor microenvironment of colorectal cancer using high-definition spatial transcriptomics"
- **Accession:** Zenodo `10.5281/zenodo.7376897` (linked from GEO `GSE226997`)
- **Files:** 16 Visium sections + 4 Visium HD captures
- **Size:** ~24 GB
- **URL:** https://zenodo.org/record/7376897  ·  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226997

---

## 2. Single-cell reference atlases

### 2.1 PanSci mouse multi-organ
- **Name:** PanSci — pan-organ mouse sci-RNA-seq atlas
- **Paper:** Zhang Z. et al., *Nature* 2024 — "A panoramic view of cell population dynamics in mammalian aging"
- **Accession:** GEO `GSE247719`
- **Files:** 3 (`counts.mtx.gz`, `barcodes.tsv.gz`, `metadata.csv.gz`)
- **Size:** ~14 GB
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247719

### 2.2 Tabula Muris Senis
- **Name:** Tabula Muris Senis (FACS + droplet)
- **Paper:** Schaum N. et al., *Nature* 2020 — "Ageing hallmarks exhibit organ-specific temporal signatures"
- **Accession:** figshare collection 4753912
- **Files:** 2 master `.h5ad` (FACS, droplet) + per-organ files (~46 files)
- **Size:** ~8.5 GB
- **URL:** https://figshare.com/projects/Tabula_Muris_Senis/64982

### 2.3 Mouse Cell Atlas 2.0 (MCA2)
- **Name:** MCA2 microwell-seq adult mouse atlas
- **Paper:** Han X. et al., *Cell* 2022 — "Mapping the Mouse Cell Atlas by Microwell-seq" (v2 update)
- **Accession:** GEO `GSE153562`
- **Files:** 1 combined `dge_rmbatch.h5ad` + per-tissue tables (~30 files)
- **Size:** ~6.8 GB
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153562

### 2.4 EasySci mouse brain aging
- **Name:** EasySci-RNA mouse brain atlas
- **Paper:** Sziraki A. et al., *Nat Genet* 2023 — "A global view of aging and Alzheimer's pathogenesis-associated cell population dynamics and molecular signatures in human and mouse brains"
- **Accession:** GEO `GSE212606`
- **Files:** 6 (matrix, barcodes, features, metadata × 2 batches)
- **Size:** ~5.2 GB
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212606

### 2.5 Hammond / Olah microglia reference
- **Name:** Human microglia single-nucleus reference
- **Paper:** Olah M. et al., *Nat Commun* 2020 — "Single cell RNA sequencing of human microglia uncovers a subset associated with Alzheimer's disease"
- **Accession:** GEO `GSE179358`
- **Files:** 1 `.h5ad` + sample metadata (~3 files)
- **Size:** ~1.4 GB
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179358

### 2.6 Burclaff 2022 human colon
- **Name:** Adult human colon single-cell atlas
- **Paper:** Burclaff J. et al., *Cell Mol Gastroenterol Hepatol* 2022 — "A proximal-to-distal survey of healthy adult human small intestine and colon epithelium by single-cell transcriptomics"
- **Accession:** Single Cell Portal `SCP1038`
- **Files:** 1 expression matrix + metadata + cluster file (3 files)
- **Size:** ~2.0 GB
- **URL:** https://singlecell.broadinstitute.org/single_cell/study/SCP1038

---

## 3. GWAS summary statistics

### 3.1 Primary trait — Alzheimer's disease
- **Name:** Bellenguez AD GWAS (Stage 1+2 meta-analysis)
- **Paper:** Bellenguez C. et al., *Nat Genet* 2022 — "New insights into the genetic etiology of Alzheimer's disease and related dementias"
- **Accession:** GWAS Catalog `GCST90027158`
- **File:** `GCST90027158_buildGRCh38.tsv.gz` (1 file)
- **Size:** ~620 MB
- **URL:** https://www.ebi.ac.uk/gwas/studies/GCST90027158  ·  ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/

### 3.2 Secondary GWAS panel (14 traits, comparator)

| # | Trait | Author / Year | GWAS Catalog accession | Size |
|---|---|---|---|---|
| 1 | Parkinson's disease (PD) | Nalls 2019 | `GCST009325` | ~110 MB |
| 2 | Amyotrophic lateral sclerosis (ALS) | van Rheenen 2021 | `GCST90027163` | ~280 MB |
| 3 | Multiple sclerosis (MS) | IMSGC 2019 | `GCST009597` | ~85 MB |
| 4 | Frontotemporal dementia (FTD) | Pottier 2019 | `GCST009666` | ~45 MB |
| 5 | Schizophrenia (SCZ) | Trubetskoy 2022 | `GCST90128468` | ~540 MB |
| 6 | Type-2 diabetes (T2D) | Mahajan 2022 | `GCST90132184` | ~720 MB |
| 7 | Bone mineral density (BMD, heel) | Morris 2019 | `GCST006979` | ~470 MB |
| 8 | C-reactive protein (CRP) | Said 2022 | `GCST90029070` | ~410 MB |
| 9 | Gastroesophageal reflux disease (GERD) | An 2019 | `GCST90000514` | ~290 MB |
| 10 | Irritable bowel syndrome (IBS) | Eijsbouts 2021 | `GCST90016564` | ~310 MB |
| 11 | Major depressive disorder (MDD) | Howard 2019 | `GCST005902` | ~250 MB |
| 12 | Attention-deficit/hyperactivity disorder (ADHD) | Demontis 2023 | `GCST90269854` | ~380 MB |
| 13 | Adult height | Yengo 2022 | `GCST90245848` | ~830 MB |
| 14 | Rheumatoid arthritis (RA) | Ishigaki 2022 | `GCST90132222` | ~520 MB |

All secondary files are pulled from the GWAS Catalog FTP under
`ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/<accession>/`.
Where the publisher posted GRCh37, we lift to GRCh38 with `CrossMap` before use.

---

## 4. eQTL & functional annotation

### 4.1 GTEx v8 eQTLs — colon transverse & sigmoid
- **Name:** GTEx v8 single-tissue cis-eQTL (Colon — Transverse; Colon — Sigmoid)
- **Paper:** GTEx Consortium, *Science* 2020 — "The GTEx Consortium atlas of genetic regulatory effects across human tissues"
- **Accession:** dbGaP `phs000424.v8.p2` (released eQTL data are open)
- **Access mode:** Live API queries, *not* bulk download
  - Endpoint: `https://gtexportal.org/api/v2/association/singleTissueEqtl`
  - Parameters used: `tissueSiteDetailId=Colon_Transverse|Colon_Sigmoid`, `gencodeId=<ENSG>`, `datasetId=gtex_v8`
- **Cached pulls:** stored under `results/eqtl_cache/{Colon_Transverse,Colon_Sigmoid}/*.json`
- **Size on disk after caching:** ~120 MB
- **URL (docs):** https://gtexportal.org/api/v2/redoc

---

## 5. gsMap resource bundle (LDSC + SNP annotations)

- **Name:** gsMap pre-built resource bundle (1KG-EUR baseline LD scores, SNP→gene mapping, plink files)
- **Paper:** Song L. et al., *Nat Genet* 2024 — "Spatially resolved mapping of cells associated with human complex traits"
- **Files:** ~1,400 files (per-chromosome LD scores, baseline annotation, SNP map)
- **Size:** ~29 GB
- **URL:** https://yanglab.westlake.edu.cn/data/gsMap/gsMap_resource.tar.gz
- **Mirror:** https://github.com/JianYang-Lab/gsMap (see README "Resource bundle")
- **Local path after extraction:** `resources/gsMap_resource/`

---

## 6. Provenance & reproducibility notes

- All accessions above were resolved on **2026-05-30**.
- For each GEO accession, we record the SRA SHA256 of the downloaded archive in
  `data_manifest/checksums.sha256` (generated by `scripts/01_fetch_all.sh`).
- GWAS Catalog files: we pin the harmonised `*_buildGRCh38.tsv.gz` where the
  curator has lifted it; otherwise we lift GRCh37 → GRCh38 with `CrossMap.py
  bed hg19ToHg38.over.chain.gz`.
- GTEx and Single Cell Portal require login-protected click-through but no
  controlled-access approval; no dbGaP application is required for any analysis
  reproduced here.
- The complete fetch script is `scripts/01_fetch_all.sh`; a dry-run that
  prints expected sizes is `scripts/01_fetch_all.sh --dry-run`.

---

*Maintainer:* steorra@stanford.edu  ·  *Last updated:* 2026-05-30
