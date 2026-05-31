## Supplementary Table Mouse 4 | Final per-gene evidence matrix (Para 5 S5–S6)

Combined evidence from **GTEx eQTL + AD-SNP overlap + coloc-ABF (full summary stats) + SharePro + SMR + HEIDI** for the 8 focal AD risk genes nominated by the colon spatial analysis. Each method has its own threshold (NG convention); the **verdict** column gives the integrated call, and **missing_analyses** is the explicit gap that the paragraph S6 calls out for CD2AP.

### 1 | Evidence matrix

| Gene | colon n_eQTL | colon best_p | AD-SNP=colon eQTL | coloc PP.H4 | coloc PP.H3 | SharePro share | SMR best p | HEIDI p | coloc✓ | SharePro✓ | SMR✓ | HEIDI✓ | Missing | Verdict |
|---|---:|---:|---|---:|---:|---:|---:|---:|:---:|:---:|:---:|:---:|---|---|
| **CR1** | 24 | 4.60e-07 | Yes | 0.9925 | 0.3416 | 0.9991 | 2.07e-06 | 0.00e+00 | ✓ | ✓ | ✓ | — | (complete) | ✓✓ shared-causal (coloc + SharePro) |
| **PICALM** | 93 | 2.49e-09 | No | 0.0001 | 0.9995 | 0.9977 | 1.34e-03 | 0.00e+00 | — | ✓ | ✓ | — | (complete) | ✓ shared-causal (multi-causal SharePro rescue) |
| **CD2AP** | 200 | 1.56e-06 | Yes | nan | nan | nan |  |  | — | — | — | — | SharePro, coloc-ABF, SMR/HEIDI | REQUIRES FOLLOW-UP: AD GWAS SNP overlaps a colon eQTL but multi-causal SharePro not yet run (SharePro + coloc-ABF + SMR/HEIDI missing) |
| **ADAM10** | 345 | 8.73e-10 | No | 0.0004 | 0.9996 | 0.0005 | 4.00e-04 | 0.00e+00 | — | — | ✓ | — | (complete) | colon eQTL present but no coloc / SharePro / SMR support |
| **APP** | 28 | 5.99e-07 | No | 0.0982 | 0.9621 | nan | 1.48e-02 | 1.47e-129 | — | — | — | — | SharePro | REQUIRES FOLLOW-UP: has colon eQTL but missing SharePro |
| **SORL1** | 0 |  | No | 0.0514 | 0.1161 | nan | 4.22e-01 | 2.61e-43 | — | — | — | — | SharePro | no colon eQTL |
| **PSEN1** | 0 |  | No | 0.037 | 0.0388 | nan | 1.49e-02 | 4.28e-17 | — | — | — | — | SharePro | no colon eQTL |
| **BIN1** | 192 | 3.20e-22 | No | 0.0 | 1.0 | nan |  |  | — | — | — | — | SharePro, SMR/HEIDI | REQUIRES FOLLOW-UP: has colon eQTL but missing SharePro + SMR/HEIDI |

### 2 | Per-method raw values (long format)

| Gene | Method | Metric | Value | Tissue |
|---|---|---|---:|---|
| CR1 | GTEx | colon_n_eQTL | 24 | Colon (Trans+Sigmoid) |
| CR1 | GTEx | colon_best_eQTL_p | 4.60e-07 | Colon (Trans+Sigmoid) |
| CR1 | GTEx | AD_SNP_is_colon_eQTL | Yes | Colon (Trans+Sigmoid) |
| CR1 | coloc-ABF | best_PP_H4 | 0.9925 | Colon (Trans+Sigmoid) |
| CR1 | coloc-ABF | best_PP_H3 | 0.3416 | Colon (Trans+Sigmoid) |
| CR1 | SharePro | share_probability | 0.9991 | Colon (multi-causal) |
| CR1 | SMR | best_SMR_p | 2.07e-06 | Colon (Trans+Sigmoid) |
| CR1 | HEIDI | p_at_best_SMR | 0.00e+00 | Colon (Trans+Sigmoid) |
| PICALM | GTEx | colon_n_eQTL | 93 | Colon (Trans+Sigmoid) |
| PICALM | GTEx | colon_best_eQTL_p | 2.49e-09 | Colon (Trans+Sigmoid) |
| PICALM | GTEx | AD_SNP_is_colon_eQTL | No | Colon (Trans+Sigmoid) |
| PICALM | coloc-ABF | best_PP_H4 | 0.0001 | Colon (Trans+Sigmoid) |
| PICALM | coloc-ABF | best_PP_H3 | 0.9995 | Colon (Trans+Sigmoid) |
| PICALM | SharePro | share_probability | 0.9977 | Colon (multi-causal) |
| PICALM | SMR | best_SMR_p | 1.34e-03 | Colon (Trans+Sigmoid) |
| PICALM | HEIDI | p_at_best_SMR | 0.00e+00 | Colon (Trans+Sigmoid) |
| CD2AP | GTEx | colon_n_eQTL | 200 | Colon (Trans+Sigmoid) |
| CD2AP | GTEx | colon_best_eQTL_p | 1.56e-06 | Colon (Trans+Sigmoid) |
| CD2AP | GTEx | AD_SNP_is_colon_eQTL | Yes | Colon (Trans+Sigmoid) |
| CD2AP | coloc-ABF | best_PP_H4 | nan | Colon (Trans+Sigmoid) |
| CD2AP | coloc-ABF | best_PP_H3 | nan | Colon (Trans+Sigmoid) |
| CD2AP | SharePro | share_probability | nan | Colon (multi-causal) |
| CD2AP | SMR | best_SMR_p |  | Colon (Trans+Sigmoid) |
| CD2AP | HEIDI | p_at_best_SMR |  | Colon (Trans+Sigmoid) |
| ADAM10 | GTEx | colon_n_eQTL | 345 | Colon (Trans+Sigmoid) |
| ADAM10 | GTEx | colon_best_eQTL_p | 8.73e-10 | Colon (Trans+Sigmoid) |
| ADAM10 | GTEx | AD_SNP_is_colon_eQTL | No | Colon (Trans+Sigmoid) |
| ADAM10 | coloc-ABF | best_PP_H4 | 0.0004 | Colon (Trans+Sigmoid) |
| ADAM10 | coloc-ABF | best_PP_H3 | 0.9996 | Colon (Trans+Sigmoid) |
| ADAM10 | SharePro | share_probability | 0.0005 | Colon (multi-causal) |
| ADAM10 | SMR | best_SMR_p | 4.00e-04 | Colon (Trans+Sigmoid) |
| ADAM10 | HEIDI | p_at_best_SMR | 0.00e+00 | Colon (Trans+Sigmoid) |
| APP | GTEx | colon_n_eQTL | 28 | Colon (Trans+Sigmoid) |
| APP | GTEx | colon_best_eQTL_p | 5.99e-07 | Colon (Trans+Sigmoid) |
| APP | GTEx | AD_SNP_is_colon_eQTL | No | Colon (Trans+Sigmoid) |
| APP | coloc-ABF | best_PP_H4 | 0.0982 | Colon (Trans+Sigmoid) |
| APP | coloc-ABF | best_PP_H3 | 0.9621 | Colon (Trans+Sigmoid) |
| APP | SharePro | share_probability | nan | Colon (multi-causal) |
| APP | SMR | best_SMR_p | 1.48e-02 | Colon (Trans+Sigmoid) |
| APP | HEIDI | p_at_best_SMR | 1.47e-129 | Colon (Trans+Sigmoid) |
| SORL1 | GTEx | colon_n_eQTL | 0 | Colon (Trans+Sigmoid) |
| SORL1 | GTEx | colon_best_eQTL_p |  | Colon (Trans+Sigmoid) |
| SORL1 | GTEx | AD_SNP_is_colon_eQTL | No | Colon (Trans+Sigmoid) |
| SORL1 | coloc-ABF | best_PP_H4 | 0.0514 | Colon (Trans+Sigmoid) |
| SORL1 | coloc-ABF | best_PP_H3 | 0.1161 | Colon (Trans+Sigmoid) |
| SORL1 | SharePro | share_probability | nan | Colon (multi-causal) |
| SORL1 | SMR | best_SMR_p | 4.22e-01 | Colon (Trans+Sigmoid) |
| SORL1 | HEIDI | p_at_best_SMR | 2.61e-43 | Colon (Trans+Sigmoid) |
| PSEN1 | GTEx | colon_n_eQTL | 0 | Colon (Trans+Sigmoid) |
| PSEN1 | GTEx | colon_best_eQTL_p |  | Colon (Trans+Sigmoid) |
| PSEN1 | GTEx | AD_SNP_is_colon_eQTL | No | Colon (Trans+Sigmoid) |
| PSEN1 | coloc-ABF | best_PP_H4 | 0.037 | Colon (Trans+Sigmoid) |
| PSEN1 | coloc-ABF | best_PP_H3 | 0.0388 | Colon (Trans+Sigmoid) |
| PSEN1 | SharePro | share_probability | nan | Colon (multi-causal) |
| PSEN1 | SMR | best_SMR_p | 1.49e-02 | Colon (Trans+Sigmoid) |
| PSEN1 | HEIDI | p_at_best_SMR | 4.28e-17 | Colon (Trans+Sigmoid) |
| BIN1 | GTEx | colon_n_eQTL | 192 | Colon (Trans+Sigmoid) |
| BIN1 | GTEx | colon_best_eQTL_p | 3.20e-22 | Colon (Trans+Sigmoid) |
| BIN1 | GTEx | AD_SNP_is_colon_eQTL | No | Colon (Trans+Sigmoid) |
| BIN1 | coloc-ABF | best_PP_H4 | 0.0 | Colon (Trans+Sigmoid) |
| BIN1 | coloc-ABF | best_PP_H3 | 1.0 | Colon (Trans+Sigmoid) |
| BIN1 | SharePro | share_probability | nan | Colon (multi-causal) |
| BIN1 | SMR | best_SMR_p |  | Colon (Trans+Sigmoid) |
| BIN1 | HEIDI | p_at_best_SMR |  | Colon (Trans+Sigmoid) |

## Notes

- **Thresholds.** coloc-ABF PP.H4 ≥ 0.8 = strong shared-causal; SharePro share ≥ 0.5 = multi-causal shared signal; SMR p < 0.0042 (Bonferroni on 12 tests); HEIDI p > 0.05 = pass (single-variant model consistent).
- **Para 5 S5 verification.** PICALM SharePro share = 0.9977 (matches 0.998 in paragraph). CR1 SharePro share = 0.9991 (matches 0.999 in paragraph). ADAM10 SharePro share = 0.0005 → multi-causal NOT supported (matches paragraph).
- **Para 5 S6 verification.** CD2AP has GTEx colon n_eQTL = 200, AD_SNP_is_colon_eQTL = Yes, but SharePro analysis was not run for CD2AP in this work (file CD2AP_v3_result.sharepro.txt does not exist). Verdict: 'REQUIRES FOLLOW-UP'. This is the gene called out by the paragraph as needing additional multi-causal colocalization work.
- **Why ADAM10 fails SharePro despite high SMR.** ADAM10 has SMR p = 4.0e-04 (Bonferroni-significant) but HEIDI p ≈ 0 (rejects single-variant model). SharePro share = 0.0005 (no shared causal under multi-causal model). Interpretation: ADAM10 has eQTL + GWAS signal but the variants are NOT shared; SMR captures a linkage artifact that HEIDI correctly flags.
- **Final causal-candidate ranking (Colon).** GOLD (coloc + SMR + HEIDI all pass): none. SILVER (coloc OR SharePro shared-causal): CR1 (both methods), PICALM (SharePro rescue). INCOMPLETE: CD2AP (SharePro not yet run). DROP: ADAM10 (multi-causal contradicts shared-causal interpretation). NO COLON EQTL: SORL1, PSEN1 (zero significant colon variants in GTEx).
- **How to extend.** To run SharePro for CD2AP: place GWAS + eQTL summary statistics + LD matrix in results/sharepro/CD2AP_v3_{gwas,eqtl,ld}.txt, then run the SharePro CLI as in `notebooks/21_eqtl_coloc_analysis.ipynb`. Output will be results/sharepro/CD2AP_v3_result.sharepro.txt and this builder will pick it up automatically on re-run.

## Summary statistics

- Focal genes: 8
- Complete pipeline (no missing analyses): 3
- Genes flagged for follow-up: 5 (CD2AP, APP, SORL1, PSEN1, BIN1)
- Genes passing coloc-ABF (PP.H4 ≥ 0.8): 1 (CR1)
- Genes passing SharePro (share ≥ 0.5): 2 (CR1, PICALM)
- Genes passing SMR (Bonferroni): 3
- Genes passing HEIDI: 0

## Data sources

| Data source | Citation | Accession | URL | Version / snapshot | Path in repo | Notes |
|---|---|---|---|---|---|---|
| **GTEx Portal v8 eQTL query** | GTEx Consortium Science 2020 | `v8` | https://gtexportal.org/home/ | queried 2026-04 (snapshot) | `results/gtex_eqtl_query_results.csv` | Per (gene × tissue) significant cis-eQTL counts + best p + AD-SNP-is-eQTL flag |
| **Bellenguez 2022 AD GWAS** | Bellenguez et al. Nat Genet 2022 | `GCST90027158` | https://www.ebi.ac.uk/gwas/studies/GCST90027158 | 2022-04 (GRCh38) | `data/gwas/Bellenguez2022_AD_withN.tsv.gz (gitignored)` | AD GWAS lead SNPs per locus |
| **coloc-ABF (full summary stats)** | Giambartolomei et al. PLoS Genet 2014 | `-` | https://github.com/chr1swallace/coloc | coloc v5.x | `results/coloc_FULL_eqtl_results.csv` | Single-causal coloc using full eQTL summary stats; PP.H0–H4 per (gene × colon tissue) |
| **SharePro multi-causal coloc** | Zhang et al. Genome Biol 2023 | `-` | https://github.com/zhwm/SharePro_coloc | SharePro v5.0.0 | `results/sharepro/{GENE}_v3_result.sharepro.txt` | Multi-causal shared-component coloc; produces "share" probability per credible set |
| **SMR + HEIDI** | Zhu et al. Nat Genet 2016 | `-` | https://yanglab.westlake.edu.cn/software/smr/ | SMR v1.3.1 | `results/smr_colon_AD_results.csv` | Summary-data-based Mendelian Randomization with HEIDI heterogeneity test |
| **Pipeline notebook** | This work | `-` | - | - | `notebooks/21_eqtl_coloc_analysis.ipynb (+ _executed)` | Orchestrates GTEx query + coloc-ABF + SharePro + SMR for the 8 focal genes |
| **Comparison data: sig-only vs full** | Supplementary Fig (Para 5 S4) | `-` | - | - | `results/coloc_sig_vs_full_comparison.csv` | Sister CSV showing PP.H4 drop when switching from sig-only to full summary |
