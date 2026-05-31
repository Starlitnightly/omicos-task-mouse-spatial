## Supplementary Table Mouse 1 | Source datasets and per-(organ × age) reference assembly

Multi-age single-cell reference atlas assembled to provide cell-type and age context for the whole-body spatial AD-risk analysis. **809,284 cells** across **126 organ–age combinations**, projected onto the 6-week whole-body Array-seq spatial reference (GSE248904) using CMap. Each organ–age sample feeds gsMap spatial-LDSC to produce per-spot AD heritability enrichment.

### 1 | Source datasets (5)

| # | Source | Paper | Accession | Assay | Organs covered | Ages offered | Preprocessing status |
|---|---|---|---|---|---|---|---|
| 1 | **PanSci** | Zhang et al. Nature 2024 | `GSE247719` | sci-RNA-seq3 (PanSci) | Kidney, Lung, Heart, Liver, Muscle, Stomach, BAT, iWAT, gWAT, Ileum, Colon, Pancreas, Skin, Thymus | 3, 6, 12, 16, 23 months | Curated cellType + Lineage + Sub_cell_type annotations from source |
| 2 | **EasySci** | Sziraki et al. Nature Aging 2023 | `GSE212606` | EasySci (sci-RNA-seq3 variant, brain-optimised) | Brain (whole) | 3, 6, 21 months | Main_cell_type + Sub_cluster from source |
| 3 | **TMS** | Schaum et al. Nature 2020 (Tabula Muris Senis) | `figshare:8273102 (TMS)` | Smart-seq2 (FACS) + 10x Chromium (Droplet) | 23 organs (BAT, Bone Marrow, Brain Myeloid/NonMyeloid, Heart, Kidney, Liver, Lung, Pancreas, Spleen, Skin, Muscle, Stomach, Colon, Thymus, Trachea, …) | 1, 3, 18, 21, 24, 30 months | cell_ontology_class from source (Tabula Muris cell ontology) |
| 4 | **MCA2** | Han et al. Cell 2022 (Mouse Cell Atlas 2.0) | `GSE153562` | Microwell-seq (MCA2) | Stomach, Ileum, Spleen, Thymus (mid-late life gaps for TMS organs) | 6, 12, 18, 21, 24 months (per-organ subset) | leiden clusters with marker-based manual annotation |
| 5 | **Hammond_microglia** | Olah et al. Nature Communications 2020 | `GSE179358` | 10x Chromium (microglia FACS-sorted, brain) | Brain (microglia only) | 6, 12, 18 months | cellType=microglial cell (single-cell-type focus) |

### 2 | Per-source contribution to the assembled reference

| Source | Samples | Distinct organs | Total cells | % of atlas | Age range |
|---|---:|---:|---:|---:|---|
| PanSci | 47 | 10 | 454,620 | 56.2% | 3–23 mo |
| TMS | 65 | 13 | 207,717 | 25.7% | 1–30 mo |
| MCA2 | 12 | 6 | 110,309 | 13.6% | 6–24 mo |
| EasySci | 3 | 1 | 26,640 | 3.3% | 3–21 mo |
| Hammond_microglia | 4 | 1 | 9,998 | 1.2% | 6–24 mo |
| **Total** | **126** | **15** | **809,284** | **100.0%** | **1–30 mo** |

### 3 | Full per-(organ × age × source) ledger (131 rows)

Multi-source rows appear where >1 dataset contributes cells to the same organ × age combination (5 such Brain samples). Each row corresponds to one CMap-projection unit feeding gsMap.

| Organ | Age (mo) | Source | Sample ID | n_cells | Source accession | Preprocessing | Mapping target |
|---|---:|---|---|---:|---|---|---|
| BAT | 3 | PanSci | `BAT_03_months` | 10,000 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Brown_Fat_CTRL1 (GSE248904) |
| BAT | 3 | TMS | `BAT_3m` | 713 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Brown_Fat_CTRL1 (GSE248904) |
| BAT | 6 | PanSci | `BAT_06_months` | 10,000 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Brown_Fat_CTRL1 (GSE248904) |
| BAT | 12 | PanSci | `BAT_12_months` | 10,000 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Brown_Fat_CTRL1 (GSE248904) |
| BAT | 16 | PanSci | `BAT_16_months` | 10,000 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Brown_Fat_CTRL1 (GSE248904) |
| BAT | 18 | TMS | `BAT_18m` | 662 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Brown_Fat_CTRL1 (GSE248904) |
| BAT | 21 | TMS | `BAT_21m` | 794 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Brown_Fat_CTRL1 (GSE248904) |
| BAT | 23 | PanSci | `BAT_23_months` | 10,000 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Brown_Fat_CTRL1 (GSE248904) |
| BAT | 24 | TMS | `BAT_24m` | 848 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Brown_Fat_CTRL1 (GSE248904) |
| BAT | 30 | TMS | `BAT_30m` | 2,015 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Brown_Fat_CTRL1 (GSE248904) |
| Bone_Marrow | 1 | TMS | `Bone_Marrow_1m` | 3,027 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Bone_Marrow_CTRL1 (GSE248904) |
| Bone_Marrow | 3 | TMS | `Bone_Marrow_3m` | 3,490 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Bone_Marrow_CTRL1 (GSE248904) |
| Bone_Marrow | 18 | TMS | `Bone_Marrow_18m` | 6,712 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Bone_Marrow_CTRL1 (GSE248904) |
| Bone_Marrow | 21 | TMS | `Bone_Marrow_21m` | 5,216 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Bone_Marrow_CTRL1 (GSE248904) |
| Bone_Marrow | 24 | TMS | `Bone_Marrow_24m` | 8,279 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Bone_Marrow_CTRL1 (GSE248904) |
| Bone_Marrow | 30 | TMS | `Bone_Marrow_30m` | 10,000 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Bone_Marrow_CTRL1 (GSE248904) |
| Brain | 3 | EasySci | `Brain_03_months` | 8,875 | `GSE212606` | Main_cell_type + Sub_cluster from source | Brain_CTRL1 (GSE248904) |
| Brain | 3 | PanSci | `Brain_3m` | 4,405 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Brain_CTRL1 (GSE248904) |
| Brain | 3 | TMS | `Brain_3m` | 4,643 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Brain_CTRL1 (GSE248904) |
| Brain | 6 | EasySci | `Brain_06_months` | 8,892 | `GSE212606` | Main_cell_type + Sub_cluster from source | Brain_CTRL1 (GSE248904) |
| Brain | 6 | Hammond_microglia | `Brain_6m` | 2,528 | `GSE179358` | cellType=microglial cell (single-cell-type focus) | Brain_CTRL1 (GSE248904) |
| Brain | 6 | PanSci | `Brain_6m` | 4,478 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Brain_CTRL1 (GSE248904) |
| Brain | 12 | Hammond_microglia | `Brain_12m` | 2,197 | `GSE179358` | cellType=microglial cell (single-cell-type focus) | Brain_CTRL1 (GSE248904) |
| Brain | 12 | MCA2 | `Brain_12m` | 4,242 | `GSE153562` | leiden clusters with marker-based manual annotation | Brain_CTRL1 (GSE248904) |
| Brain | 18 | Hammond_microglia | `Brain_18m` | 2,614 | `GSE179358` | cellType=microglial cell (single-cell-type focus) | Brain_CTRL1 (GSE248904) |
| Brain | 18 | TMS | `Brain_18m` | 4,022 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Brain_CTRL1 (GSE248904) |
| Brain | 21 | EasySci | `Brain_21_months` | 8,873 | `GSE212606` | Main_cell_type + Sub_cluster from source | Brain_CTRL1 (GSE248904) |
| Brain | 21 | TMS | `Brain_21m` | 11,824 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Brain_CTRL1 (GSE248904) |
| Brain | 24 | Hammond_microglia | `Brain_24m` | 2,659 | `GSE179358` | cellType=microglial cell (single-cell-type focus) | Brain_CTRL1 (GSE248904) |
| Brain | 24 | TMS | `Brain_24m` | 2,677 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Brain_CTRL1 (GSE248904) |
| Colon | 3 | PanSci | `Colon_03_months` | 10,000 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Colon_CTRL1 (GSE248904) |
| Colon | 3 | TMS | `Colon_3m` | 3,987 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Colon_CTRL1 (GSE248904) |
| Colon | 6 | PanSci | `Colon_06_months` | 10,000 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Colon_CTRL1 (GSE248904) |
| Colon | 12 | PanSci | `Colon_12_months` | 10,000 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Colon_CTRL1 (GSE248904) |
| Colon | 16 | PanSci | `Colon_16_months` | 10,000 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Colon_CTRL1 (GSE248904) |
| Colon | 18 | TMS | `Colon_18m` | 2,369 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Colon_CTRL1 (GSE248904) |
| Colon | 21 | MCA2 | `Colon_21m` | 10,000 | `GSE153562` | leiden clusters with marker-based manual annotation | Colon_CTRL1 (GSE248904) |
| Colon | 23 | PanSci | `Colon_23_months` | 10,000 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Colon_CTRL1 (GSE248904) |
| Colon | 24 | TMS | `Colon_24m` | 1,955 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Colon_CTRL1 (GSE248904) |
| Colon | 30 | TMS | `Colon_30m` | 1,887 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Colon_CTRL1 (GSE248904) |
| Heart | 1 | TMS | `Heart_1m` | 946 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Heart_CTRL1 (GSE248904) |
| Heart | 3 | PanSci | `Heart_03_months` | 10,000 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Heart_CTRL1 (GSE248904) |
| Heart | 3 | TMS | `Heart_3m` | 534 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Heart_CTRL1 (GSE248904) |
| Heart | 6 | PanSci | `Heart_06_months` | 10,000 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Heart_CTRL1 (GSE248904) |
| Heart | 12 | PanSci | `Heart_12_months` | 10,000 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Heart_CTRL1 (GSE248904) |
| Heart | 16 | PanSci | `Heart_16_months` | 10,000 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Heart_CTRL1 (GSE248904) |
| Heart | 18 | TMS | `Heart_18m` | 1,143 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Heart_CTRL1 (GSE248904) |
| Heart | 21 | TMS | `Heart_21m` | 918 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Heart_CTRL1 (GSE248904) |
| Heart | 23 | PanSci | `Heart_23_months` | 10,000 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Heart_CTRL1 (GSE248904) |
| Heart | 24 | TMS | `Heart_24m` | 768 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Heart_CTRL1 (GSE248904) |
| Heart | 30 | TMS | `Heart_30m` | 4,304 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Heart_CTRL1 (GSE248904) |
| Ileum | 3 | PanSci | `Ileum_03_months` | 9,994 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Small_Intestine_CTRL1 (GSE248904) |
| Ileum | 6 | PanSci | `Ileum_06_months` | 9,994 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Small_Intestine_CTRL1 (GSE248904) |
| Ileum | 12 | MCA2 | `Ileum_12m` | 9,981 | `GSE153562` | leiden clusters with marker-based manual annotation | Small_Intestine_CTRL1 (GSE248904) |
| Ileum | 12 | PanSci | `Ileum_12_months` | 9,996 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Small_Intestine_CTRL1 (GSE248904) |
| Ileum | 16 | PanSci | `Ileum_16_months` | 9,995 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Small_Intestine_CTRL1 (GSE248904) |
| Ileum | 18 | MCA2 | `Ileum_18m` | 9,983 | `GSE153562` | leiden clusters with marker-based manual annotation | Small_Intestine_CTRL1 (GSE248904) |
| Ileum | 21 | MCA2 | `Ileum_21m` | 9,985 | `GSE153562` | leiden clusters with marker-based manual annotation | Small_Intestine_CTRL1 (GSE248904) |
| Ileum | 23 | PanSci | `Ileum_23_months` | 9,998 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Small_Intestine_CTRL1 (GSE248904) |
| Ileum | 24 | MCA2 | `Ileum_24m` | 9,986 | `GSE153562` | leiden clusters with marker-based manual annotation | Small_Intestine_CTRL1 (GSE248904) |
| Kidney | 1 | TMS | `Kidney_1m` | 2,422 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Kidney_CTRL1 (GSE248904) |
| Kidney | 3 | PanSci | `Kidney_03_months` | 9,802 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Kidney_CTRL1 (GSE248904) |
| Kidney | 3 | TMS | `Kidney_3m` | 2,254 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Kidney_CTRL1 (GSE248904) |
| Kidney | 6 | PanSci | `Kidney_06_months` | 9,773 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Kidney_CTRL1 (GSE248904) |
| Kidney | 12 | PanSci | `Kidney_12_months` | 9,795 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Kidney_CTRL1 (GSE248904) |
| Kidney | 16 | PanSci | `Kidney_16_months` | 9,800 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Kidney_CTRL1 (GSE248904) |
| Kidney | 18 | TMS | `Kidney_18m` | 3,040 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Kidney_CTRL1 (GSE248904) |
| Kidney | 21 | TMS | `Kidney_21m` | 2,204 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Kidney_CTRL1 (GSE248904) |
| Kidney | 23 | PanSci | `Kidney_23_months` | 9,795 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Kidney_CTRL1 (GSE248904) |
| Kidney | 24 | TMS | `Kidney_24m` | 5,508 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Kidney_CTRL1 (GSE248904) |
| Kidney | 30 | TMS | `Kidney_30m` | 5,511 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Kidney_CTRL1 (GSE248904) |
| Liver | 1 | TMS | `Liver_1m` | 2,785 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Liver_CTRL1 (GSE248904) |
| Liver | 3 | PanSci | `Liver_03_months` | 9,981 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Liver_CTRL1 (GSE248904) |
| Liver | 3 | TMS | `Liver_3m` | 729 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Liver_CTRL1 (GSE248904) |
| Liver | 6 | PanSci | `Liver_06_months` | 9,973 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Liver_CTRL1 (GSE248904) |
| Liver | 12 | PanSci | `Liver_12_months` | 9,969 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Liver_CTRL1 (GSE248904) |
| Liver | 16 | PanSci | `Liver_16_months` | 9,969 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Liver_CTRL1 (GSE248904) |
| Liver | 18 | TMS | `Liver_18m` | 1,176 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Liver_CTRL1 (GSE248904) |
| Liver | 21 | TMS | `Liver_21m` | 367 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Liver_CTRL1 (GSE248904) |
| Liver | 23 | PanSci | `Liver_23_months` | 9,983 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Liver_CTRL1 (GSE248904) |
| Liver | 24 | TMS | `Liver_24m` | 945 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Liver_CTRL1 (GSE248904) |
| Liver | 30 | TMS | `Liver_30m` | 2,272 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Liver_CTRL1 (GSE248904) |
| Lung | 1 | TMS | `Lung_1m` | 2,507 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Lung_CTRL1 (GSE248904) |
| Lung | 3 | PanSci | `Lung_03_months` | 9,996 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Lung_CTRL1 (GSE248904) |
| Lung | 3 | TMS | `Lung_3m` | 3,610 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Lung_CTRL1 (GSE248904) |
| Lung | 6 | PanSci | `Lung_06_months` | 9,991 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Lung_CTRL1 (GSE248904) |
| Lung | 12 | PanSci | `Lung_12_months` | 9,995 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Lung_CTRL1 (GSE248904) |
| Lung | 16 | PanSci | `Lung_16_months` | 9,996 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Lung_CTRL1 (GSE248904) |
| Lung | 18 | TMS | `Lung_18m` | 4,376 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Lung_CTRL1 (GSE248904) |
| Lung | 21 | TMS | `Lung_21m` | 2,950 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Lung_CTRL1 (GSE248904) |
| Lung | 23 | PanSci | `Lung_23_months` | 9,997 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Lung_CTRL1 (GSE248904) |
| Lung | 24 | TMS | `Lung_24m` | 1,870 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Lung_CTRL1 (GSE248904) |
| Lung | 30 | TMS | `Lung_30m` | 9,994 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Lung_CTRL1 (GSE248904) |
| Muscle | 3 | PanSci | `Muscle_03_months` | 9,995 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Muscle_CTRL1 (GSE248904) |
| Muscle | 3 | TMS | `Muscle_3m` | 1,102 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Muscle_CTRL1 (GSE248904) |
| Muscle | 6 | PanSci | `Muscle_06_months` | 9,999 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Muscle_CTRL1 (GSE248904) |
| Muscle | 12 | PanSci | `Muscle_12_months` | 9,995 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Muscle_CTRL1 (GSE248904) |
| Muscle | 16 | PanSci | `Muscle_16_months` | 9,997 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Muscle_CTRL1 (GSE248904) |
| Muscle | 18 | TMS | `Muscle_18m` | 1,521 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Muscle_CTRL1 (GSE248904) |
| Muscle | 21 | TMS | `Muscle_21m` | 3,528 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Muscle_CTRL1 (GSE248904) |
| Muscle | 23 | PanSci | `Muscle_23_months` | 9,995 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Muscle_CTRL1 (GSE248904) |
| Muscle | 24 | TMS | `Muscle_24m` | 1,231 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Muscle_CTRL1 (GSE248904) |
| Muscle | 30 | TMS | `Muscle_30m` | 3,284 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Muscle_CTRL1 (GSE248904) |
| Pancreas | 3 | TMS | `Pancreas_3m` | 1,588 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Pancreas_CTRL1 (GSE248904) |
| Pancreas | 18 | TMS | `Pancreas_18m` | 2,071 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Pancreas_CTRL1 (GSE248904) |
| Pancreas | 21 | TMS | `Pancreas_21m` | 1,437 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Pancreas_CTRL1 (GSE248904) |
| Pancreas | 30 | TMS | `Pancreas_30m` | 2,642 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Pancreas_CTRL1 (GSE248904) |
| Skin | 3 | TMS | `Skin_3m` | 2,346 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Skin_CTRL1 (GSE248904) |
| Skin | 21 | TMS | `Skin_21m` | 4,352 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Skin_CTRL1 (GSE248904) |
| Skin | 24 | TMS | `Skin_24m` | 1,122 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Skin_CTRL1 (GSE248904) |
| Spleen | 1 | TMS | `Spleen_1m` | 2,801 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Spleen_CTRL1 (GSE248904) |
| Spleen | 3 | TMS | `Spleen_3m` | 6,480 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Spleen_CTRL1 (GSE248904) |
| Spleen | 6 | MCA2 | `Spleen_6m` | 9,412 | `GSE153562` | leiden clusters with marker-based manual annotation | Spleen_CTRL1 (GSE248904) |
| Spleen | 12 | MCA2 | `Spleen_12m` | 9,569 | `GSE153562` | leiden clusters with marker-based manual annotation | Spleen_CTRL1 (GSE248904) |
| Spleen | 18 | TMS | `Spleen_18m` | 5,929 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Spleen_CTRL1 (GSE248904) |
| Spleen | 21 | TMS | `Spleen_21m` | 6,186 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Spleen_CTRL1 (GSE248904) |
| Spleen | 24 | TMS | `Spleen_24m` | 3,954 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Spleen_CTRL1 (GSE248904) |
| Spleen | 30 | TMS | `Spleen_30m` | 8,615 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Spleen_CTRL1 (GSE248904) |
| Stomach | 3 | PanSci | `Stomach_03_months` | 9,440 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Stomach_CTRL1 (GSE248904) |
| Stomach | 6 | PanSci | `Stomach_06_months` | 9,361 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Stomach_CTRL1 (GSE248904) |
| Stomach | 12 | MCA2 | `Stomach_12m` | 8,218 | `GSE153562` | leiden clusters with marker-based manual annotation | Stomach_CTRL1 (GSE248904) |
| Stomach | 12 | PanSci | `Stomach_12_months` | 9,457 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Stomach_CTRL1 (GSE248904) |
| Stomach | 16 | PanSci | `Stomach_16_months` | 9,353 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Stomach_CTRL1 (GSE248904) |
| Stomach | 18 | MCA2 | `Stomach_18m` | 9,552 | `GSE153562` | leiden clusters with marker-based manual annotation | Stomach_CTRL1 (GSE248904) |
| Stomach | 23 | PanSci | `Stomach_23_months` | 9,353 | `GSE247719` | Curated cellType + Lineage + Sub_cell_type annotations from source | Stomach_CTRL1 (GSE248904) |
| Stomach | 24 | MCA2 | `Stomach_24m` | 9,381 | `GSE153562` | leiden clusters with marker-based manual annotation | Stomach_CTRL1 (GSE248904) |
| Thymus | 3 | TMS | `Thymus_3m` | 1,430 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Thymus_CTRL1 (GSE248904) |
| Thymus | 12 | MCA2 | `Thymus_12m` | 10,000 | `GSE153562` | leiden clusters with marker-based manual annotation | Thymus_CTRL1 (GSE248904) |
| Thymus | 18 | TMS | `Thymus_18m` | 2,040 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Thymus_CTRL1 (GSE248904) |
| Thymus | 21 | TMS | `Thymus_21m` | 2,792 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Thymus_CTRL1 (GSE248904) |
| Thymus | 24 | TMS | `Thymus_24m` | 3,013 | `figshare:8273102 (TMS)` | cell_ontology_class from source (Tabula Muris cell ontology) | Thymus_CTRL1 (GSE248904) |

## Notes

- **CMap projection.** Each single-cell sample is mapped onto the whole-body spatial reference (GSE248904 CTRL replicate 1) using CMap classification + spot assignment. See scripts/preprocess_new_sc.py and scripts/process_remaining_gaps.py. A spot-classification confidence (CMAP_quality) is retained per cell.
- **gsMap input.** Each (organ × age × source) projected h5ad is one input to `gsmap quick_mode`; outputs include per-spot AD -log10(p) and per-annotation Cauchy combination.
- **Multi-source Brain samples.** 5 Brain samples (Brain_3m, Brain_6m, Brain_12m, Brain_18m, Brain_24m) each combine ≥2 sources to maximise cell-type coverage: PanSci/EasySci supply non-microglia neurons; Hammond and TMS add microglia depth; MCA2 fills aged gaps. See Fig. 4c pie-matrix in notebooks/paper_figure.ipynb.
- **Why 15 organs and not 16.** Lymph_Node spatial spots exist in GSE248904 but no aged single-cell dataset covers lymph node directly; lymph node was excluded from the multi-age atlas (still appears in the 6-week whole-body baseline).
- **Why 1–30 mo and not 1–32 mo.** TMS droplet covers up to 30 m; no public single-cell mouse atlas reliably covers 31–32 months for the organs analysed here.
- **Provenance script.** scripts/build_age_provenance_per_source.py rebuilds results/age_sample_provenance_per_source.csv from raw obs metadata + CMap output paths.

## Summary statistics

- Total single-cell sources: 5
- Total samples (organ × age × source rows): 131
- Distinct (organ × age) combinations: 126
- Distinct organs: 15
- Distinct ages: 10 (1, 3, 6, 12, 16, 18, 21, 23, 24, 30 months)
- Total cells: 809,284
- Largest source: PanSci (454,620 cells, 56.2%)
- Smallest source: Hammond_microglia (9,998 cells)
- Multi-source samples: 5 (all Brain)

## Data sources

| Data source | Citation | Accession | URL | Version / snapshot | Path in repo | Notes |
|---|---|---|---|---|---|---|
| **PanSci** | Zhang et al. Nature 2024 | `GSE247719` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247719 | 2024-06 | `data/pansci/ (gitignored, see data_manifest/README.md §2.1)` | sci-RNA-seq3 mouse pan-organ aging atlas |
| **EasySci brain** | Sziraki et al. Nature Aging 2023 | `GSE212606` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212606 | 2023-08 | `data/pansci/brain/ (gitignored)` | EasySci, brain-optimised |
| **TMS (Tabula Muris Senis)** | Schaum et al. Nature 2020 | `figshare:8273102` | https://figshare.com/articles/dataset/8273102 | 2020-07 | `data/tms/per_organ/ (gitignored)` | Smart-seq2 FACS + 10x Chromium Droplet, 23 organs × 6 ages |
| **MCA2 (Mouse Cell Atlas 2.0)** | Han et al. Cell 2022 | `GSE153562` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153562 | 2022-09 | `data/mca2/ (gitignored)` | Microwell-seq, fills aged-organ gaps |
| **Hammond microglia** | Olah et al. Nat Commun 2020 | `GSE179358` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179358 | 2020-04 | `data/microglia_aging/ (gitignored)` | 10x microglia FACS-sorted |
| **Whole-body Array-seq mouse** | Clevenger et al. Cell 2026 | `GSE248904` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248904 | 2026-02 | `data/st/ + data/st/per_organ/ (gitignored)` | Spatial reference for CMap projection |
| **CMap (single-cell → spatial mapping)** | scripts/preprocess_new_sc.py + scripts/process_remaining_gaps.py | `-` | - | - | `models/cmap_output/{organ}/{sample}.h5ad (gitignored)` | Projects each source single-cell sample onto Clevenger 2026 spatial reference |
| **gsMap** | Liu et al. Nature Methods 2024 | `GitHub: JianYang-Lab/gsMap` | https://github.com/JianYang-Lab/gsMap | v1.x | `models/gsmap_age_output/{sample}/ (gitignored)` | Spatial S-LDSC per (sample × annotation) Cauchy combination |
| **Derived provenance CSV** | This work | `-` | - | - | `results/age_sample_provenance_per_source.csv` | Rebuilt by scripts/build_age_provenance_per_source.py |
