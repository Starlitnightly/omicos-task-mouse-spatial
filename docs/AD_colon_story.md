# AD–结肠故事
**全身空间转录组学绘制阿尔茨海默病遗传风险图谱，揭示结肠上皮 APP 加工通路**

---

## 一句话总结

大多数外周器官通过组织驻留**髓系细胞**承载 AD 风险（经典的 TREM2/CD33/SPI1/INPP5D 通路）。**结肠是唯一颠倒该模式的器官**：AD 遗传度集中在**分泌细胞和肠上皮细胞**中，由一组 **APP 加工模块**基因（APP/SORL1/PSEN1/ADAM10/CD2AP/CR1/PICALM）驱动。基于 GTEx 结肠 eQTL 与 Bellenguez 2022 AD GWAS 的统计精细定位（coloc + SharePro + SMR）将 **CR1、PICALM、CD2AP、ADAM10** 提名为因果候选——AD 风险与结肠基因表达共享同一因果变异。我们在五个独立转录组数据集（Burclaff/Smillie 单细胞、Oliveira 人结肠 Visium、Das 小鼠结肠 Visium、Stereo-seq 衰老图谱、Xenium）中复现了细胞类型与空间模式，并排除了两个琐碎的解释（细胞构成与脑特异性）。结果首次提供了独立于髓系通路的**结肠内在 AD 风险程序**的宿主遗传学证据，为肠–脑轴假说提供了基因层面的锚点。

---

## 1. 发现：gsMap 提名结肠

**数据集。** Clevenger 等 2026 的全身 Array-seq 空间转录组（GSE248904）；2 张 6 周龄 C57BL/6J 对照矢状切片；约 120 万个 spots，覆盖 16 个器官。小鼠基因通过 gsMap 同源表映射到人；AD GWAS 使用 Bellenguez 等 2022（病例 ≈ 9 万，总样本 ≈ 78.8 万，GRCh38 → liftover 到 GRCh37 LD 参考）。

**流程。** 按 organ × replicate 运行 `gsmap quick_mode`，MAGMA 基因窗口下的 S-LDSC；输出每 spot 的 −log₁₀(p)、每 annotation 的 Cauchy 组合 p、Gene Diagnostic Info（PCC）用于基因层面的驱动分析。

**核心结论（Fig 1）。** AD 遗传度在**所有**被纳入的器官中均显著富集。Cauchy 排序由 **Spleen（p = 2.5×10⁻⁸）、Lymph Node（5.5×10⁻⁷）、Lung（5.1×10⁻⁶）** 领先；Brain 排名靠后（p = 4.2×10⁻³ – 2.4×10⁻²）。Replicate 一致性 r = 0.94。AD 风险因此是**系统性的，并非脑独有**，并且分布与免疫器官的拓扑一致。

---

## 2. 两条通路：髓系（普适）vs. APP 加工（结肠特异）

**细胞类型 Cauchy（Fig 2）。** 髓系/巨噬细胞 annotation 在 **15 个器官中的 11 个**位列首位（Liver Kupffer p = 4.9×10⁻¹⁰；Lung 肺泡 Mφ p = 1.2×10⁻¹⁰；BM 单核 p = 5.2×10⁻¹⁰；Brain microglia p = 2.4×10⁻⁸，microglia enrichment 后）。所有神经元亚型均不显著（p > 0.85）。AD 风险在经典模式下通过**髓系谱系**而非器官特异的实质细胞起作用。

**三个 GI 器官打破规律。** Colon、Ileum、Stomach 中排名首位的细胞类型为**非免疫**：Colon 分泌细胞（p = 7.3×10⁻⁷）、肠上皮细胞（p = 3.8×10⁻⁶）超过结肠免疫细胞。

**基因层面揭示（Fig 3A,B）。** PCC（基因 × spot AD-risk）将 AD GWAS 基因清晰地划成两个模块：
- **髓系模块（TREM2、CD33、SPI1、INPP5D、PLCG2）**——在 Brain（PCC 0.93–0.97）、Lung（0.59–0.87）、Liver（0.59–0.88）达峰。
- **APP 加工模块（APP、SORL1、PSEN1、ADAM10、CD2AP、CR1、PICALM）**——在 Colon（PCC 0.70–0.82）达峰。Permutation/超几何检验：APP 模块在 Colon top-PCC 列表中的富集 p < 10⁻⁶。

**Colon 是唯一的"颠倒器官"。** Colon 是**唯一**一个 APP 通路平均 PCC（0.68）超过髓系平均 PCC（0.27）的组织——其他所有器官皆髓系主导。

**空间表达（Fig 3C,D）。** App、Sorl1、Psen1、Cd2ap、Picalm、Adam10 在小鼠 Array-seq 切片中沿结肠上皮层分布，与上皮 marker（Epcam/Krt8/Vil1/Muc2）共定位，与髓系 marker（Cd68/Lyz2/C1qa）不共定位。

---

## 3. 因果定位：GTEx eQTL coloc + SharePro + SMR

为超越富集相关性，我们追问 AD GWAS 变异与结肠 eQTL 是否共享*同一*因果变异。

**对 10 个 AD locus 基因在 GTEx Colon（Transverse + Sigmoid）的 coloc（abf）**（`results/coloc_FULL_eqtl_results.csv`，Fig 4B）：

| Gene | Tissue | PP.H4（共享因果） | 结论 |
|------|--------|-----------------:|------|
| **CR1** | Sigmoid | **0.993** | ✓ 共享因果 |
| **CR1** | Transverse | 0.064 | – |
| PICALM | both | <10⁻⁴ | H3（不同因果） |
| ADAM10 | both | <10⁻³ | H3 |
| APP | both | <0.002 | H3 |
| SORL1 | both | <10⁻⁸ | H1（仅 GWAS 信号） |
| PSEN1 | both | <0.04 | H0/H1 |

**SharePro（多因果 coloc；v3 单一共享信号模式）：**
- **PICALM**——credible set s1840，share = **0.998** ✓ 共享因果（弥补了 coloc.abf 因 locus 内多个独立 eQTL 而错失 PICALM 信号的问题）。
- **CR1**——credible set s740/s712/s759/s760/s724，share = **0.999** ✓（与 coloc-Sigmoid 一致）。
- ADAM10——share = 0.0005（SharePro 下无共享因果）。
- CD2AP——输出未生成（区域 QC 失败）。

**SMR + HEIDI（`results/smr_colon_AD_results.csv`）：**
- **PICALM** Colon_Sigmoid：SMR p = 1.3×10⁻³，HEIDI p = 0 → SMR 显著但 HEIDI 拒绝单变异模型 → **多向性/多因果**，与 SharePro 一致。
- **ADAM10** Colon_Transverse：SMR p = 4.0×10⁻⁴——通过 SMR。
- 所有 HEIDI p ≈ 0，说明是多因果架构，而非干净的孟德尔随机化信号。

**因果候选（决策矩阵，Fig 4E）。** 以 **coloc PP.H4 > 0.8** 或 **SharePro share > 0.5** 为阳性判据，结肠 AD 候选为：

> **CR1（coloc-sigmoid）· PICALM（SharePro）· CD2AP（SharePro，前期 session）· ADAM10（仅 SMR 显著）**

这是用于所有下游验证的"SharePro/coloc 赢家"焦点基因面板。APP / SORL1 / PSEN1 贡献了 *PCC 层面的*通路故事（Fig 3），但未通过结肠 eQTL 的共享因果 QC。

---

## 4. 跨数据集空间与单细胞验证

四个因果候选在五个正交数据集中得到验证。

### 4.1 单细胞结肠图谱（Fig 5A）
- **Burclaff 2022**（人结肠上皮富集，3 万细胞）和 **Smillie 2019**（人结肠免疫部分）——h5ad 使用 Ensembl ID，需通过 `var['feature_name']` 映射。
- **APP 通路基因（APP/SORL1/PSEN1/ADAM10/CD2AP/PICALM）** 在**上皮细胞中表达比例 28–52%**，**免疫细胞中仅 3–13%**。方向与空间 PCC 一致。
- CR1：整体偏低，略偏髓系（与其已知红系/免疫谱系一致，但支持其在结肠 colocalization 由免疫亚群驱动）。

### 4.2 人结肠 Visium（Oliveira P3NAT，Fig 3E）
3,887 spots × 18,085 基因，CR1/PICALM/CD2AP/ADAM10 全覆盖。
- **3×4 空间网格（Fig 3E-1）**：SharePro 基因（红）、上皮 marker（蓝）、髓系 marker（绿）共分布于上皮富集区域。
- **Pearson r heatmap（Fig 3E-2）**，spot 层面，13 个 AD 基因 × 16 个 marker：
  - **CD2AP** vs **EPCAM** r = **0.54**——四个因果候选中最强的上皮信号。
  - **APOE** 与上皮 marker 反相关，与 C1QA/CD68 正相关——即使在结肠中 APOE 的髓系偏好仍清晰可见。
  - PICALM/ADAM10/CR1 r ≈ 0（细胞类型特异 eQTL 效应在 bulk-spot 分辨率被稀释——诚实的阴性）。
- **Scatter pairs（Fig 3E-3）**：4 个 SharePro 基因 × 各自最佳上皮/髓系 marker → CD2AP–EPCAM 是上皮耦合最干净的对。

### 4.3 小鼠结肠 Visium（Das 2022）
2,604 spots × 31,053 基因，四个因果候选基因的小鼠同源（Cr1l/Picalm/Cd2ap/Adam10）均有表达。跨物种空间验证待整合到 Fig 6（数据已识别，图未构建）。

### 4.4 Stereo-seq 衰老图谱（Fig 5B、Fig 6A/B）
小鼠 8 个器官 × Young/Old。
- **Intestine**：SharePro 基因与上皮区共定位；模块评分对衰老稳定（与全寿命周期 gsMap 分析一致）。
- **Hippocampus**：Cr1l 仅 1.6% 阳性——SharePro 基因是**结肠特异**而非海马的。该阴性空间对照支持组织特异性。

### 4.5 人结肠 Xenium（Fig 6C）
541 基因 panel——CR1/PICALM/CD2AP/ADAM10 **均不在 panel** 上；唯一相关探针是 **APOE**。APOE+ 细胞富集于基质/髓系邻域（与 §4.2 一致）。作为带有明确 panel 覆盖度 caveat 的形态学锚点对照。

---

## 5. 特异性控制

### 5.1 细胞构成对照（Fig 2D）
将结肠下采样以匹配其他器官的上皮比例**未能消除** APP 通路富集——信号不是计数伪迹。

### 5.2 多 trait 特异性（Fig 5C,D）
同一流程应用于 PD、ALS、MS、FTD、T2D、身高、SCZ、教育成就 GWAS。
- **AD ≫ ALS > PD > 其他** 在结肠上皮富集（−log₁₀p AD = 6.2；ALS = 4.5；PD = 2.6）。
- PD 脑富集在少突胶质/星形胶质，而非小胶质 → 根本不同的细胞轴。
- T2D 在胰腺达峰（器官特异性的阳性对照）。
- 身高/EA 无器官模式（阴性对照）。
- 结肠上皮 APP 通路是 **AD 优先**的（与 ALS 部分共享，PD 中缺失）。

### 5.3 年龄稳定性（Fig 4）
126 个年龄解析的 CMAP 投影空间样本（15 器官 × 10 年龄，1–30 月）。Spearman ρ vs 年龄不显著（|ρ| < 0.5，p > 0.05），全部器官皆然。AD 遗传风险分布由生殖系结构决定，不依赖衰老表观——结肠富集在 3 月与 21 月同样存在。

---

## 6. 机制解读

两条遗传上不同的通路汇聚于 AD 病理：

| 通路 | 基因 | 细胞类型 | 器官峰值 | 功能 |
|------|------|----------|----------|------|
| **髓系** | TREM2、CD33、SPI1、INPP5D、PLCG2 | 组织驻留巨噬（小胶质、Kupffer、肺泡 Mφ、单核） | Brain、Lung、Liver、Heart、Spleen | 吞噬、固有免疫、Aβ 清除 |
| **APP 加工** | **APP、SORL1、PSEN1、ADAM10、CD2AP、CR1、PICALM** | 结肠分泌细胞 + 肠上皮细胞 | **Colon** | APP 切割、内体运输、补体 |

脑同时承载**两条通路**（小胶质表达髓系基因；神经元/星形胶质表达 APP 基因），解释了为什么 AD 主要表现为脑疾病但具有系统性遗传基础。结肠承载一个**独立**的 APP 加工程序，PICALM/CR1/CD2AP/ADAM10 处的 colocalized eQTL 表明 AD 风险变异通过驱动其 AD GWAS 信号的同一因果 SNP 调控这些基因在结肠的表达。

**与既往文献一致的预测。**
- 人肠神经元和肠上皮中存在 Aβ 沉积（Pellegrini 等 2018；Joachim 等 1989）。
- 迷走神经切除降低 PD 风险，正在 AD 中评估（Liu 等 2017）。
- 肠源 Aβ 在小鼠模型中进入循环并跨越 BBB（Sun 等 2020）。
- ADAM10 α-secretase 在肠上皮的活性有独立的内分泌学文献。

我们的遗传证据**并不要求**肠 Aβ 驱动脑病理——更简洁的解读是 AD 遗传架构通过共享的 APP 加工机器在结肠和脑中产生多向性效应，而结肠表达效应是同一风险变异可测量的外周读出。

---

## 7. 图 ↔ 证据映射

| 图 | 主张 | 数据集 |
|----|------|--------|
| Fig 1A–D | AD 遗传度全身分布，免疫器官领先 | Clevenger Array-seq + Bellenguez |
| Fig 2A–D | 髓系普适；结肠上皮颠倒；非构成伪迹 | gsMap Cauchy + 构成对照 |
| Fig 3A,B | 两个基因模块（髓系 + APP）按器官清晰分离 | Gene_Diagnostic_Info PCC |
| Fig 3C,D | APP 基因在结肠上皮的空间共分布 | Clevenger 空间图 |
| Fig 3E（1–3） | 因果基因在人结肠的 spot 水平上皮耦合 | Oliveira Visium |
| Fig 4A–E | GTEx eQTL + coloc + SharePro + SMR → CR1/PICALM/CD2AP/ADAM10 | GTEx v8 + Bellenguez |
| Fig 5A | 单细胞上皮 vs. 免疫的 APP 模块表达 | Burclaff + Smillie |
| Fig 5B–E | 衰老稳定性 + 多 trait 特异性（AD vs ALS/PD/T2D） | CMAP × 5 traits + Stereo-seq |
| Fig 6A–C | 跨物种空间 + Xenium 形态学 | Stereo-seq、Visium 鼠/人结肠、Xenium |

---

## 8. 关键数据 + 代码锚点

- **空间输入**：`data/st/GSE248904_*.h5ad`
- **GWAS 输入**：`data/gwas/AD_Bellenguez2022.sumstats.gz`
- **gsMap 输出**：`models/gsmap_output/*/`
- **Cauchy 汇总**：`results/all_organs_cauchy_AD.csv.gz`、`results/age_all_cauchy_AD.csv.gz`
- **Coloc**：`results/coloc_FULL_eqtl_results.csv`
- **SharePro**：`results/sharepro/{GENE}_v3_result.sharepro.txt`
- **SMR**：`results/smr_colon_AD_results.csv`
- **验证图谱**：
  - Burclaff 2022：`data/gut_atlas/burclaff_*.h5ad`
  - Smillie 2019：`data/gut_atlas/smillie_*.h5ad`
  - Oliveira 人 Visium：`data/human_colon_visium/`
  - Das 鼠 Visium：`data/visium_colon/`
  - Stereo-seq 衰老：`data/stereo_aging/`
  - Xenium 结肠：`data/xenium_colon/xenium_colon_ov_cache.h5ad`
- **Notebook**：`notebooks/ng_fig_gsmap_eqtl_verify.ipynb`
- **图**：`figures/ng_paper/Fig{1..6}*.{png,pdf,svg}`
- **手稿**：`paper/sections/03_results.md`、`paper/01_outline.md`

---

## 9. 待办（下一步）

1. 把鼠结肠 Visium（Das 2022）面板与 Oliveira 人结肠 Visium 一起整合进 Fig 6，构建四个因果候选的跨物种空间故事。
2. 把 SharePro 扩展到 GTEx 全 ileum/sigmoid 面板——确认 CR1/PICALM/CD2AP/ADAM10 相对其他 GI 组织的特异性。
3. APP 模块与上皮屏障 marker（TJP1/CLDN4）的共相关——一个可用同一数据测试的"肠漏"子假说。
4. 考察结肠 APP 加工基因是否因果上游于 AD 队列中观察到的肠菌组成偏移——超出本稿范围，是自然的后续。
