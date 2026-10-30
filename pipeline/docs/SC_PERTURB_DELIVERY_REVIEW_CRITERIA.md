# scPerturb 客户交付审核标准

版本：1.1；生效日期：2026-10-20。

## 判定规则

全部必审项通过为 PASS；任一未达标为 NOT_PASS；无未达标但必审信息不足为 PENDING。后两者不能标记为验收合格。分析结果继续保留。
使用源计数或原生最高可用精度判定，边界包含在内；缺失、零分母和非有限值不算通过。改变门槛须更新版本并保留历史审核，不改变分析数据或统计定义。

## 定量 QC 门槛

| 编号 | 类别 | 指标 | 通过条件 | 统计口径 |
|---|---|---|---|---|
| QC01 | 测序深度 | GEX input reads / passing cells | >= 10,000.00 | GEX reads / passing cells |
| QC02 | 测序深度 | CRISPR input reads / passing cells | >= 800.00 | CRISPR reads / passing cells |
| QC03 | 碱基质量 | 总体 Q30 bases 比例 | >= 80.00% | GEX+CRISPR 的 Q30 bases / 总 bases |
| QC04 | GEX | Barcode 匹配率 | >= 80.00% | Total barcoded reads / Total input reads |
| QC05 | GEX | Transcriptome mapping rate | >= 70.00% | Mapped reads / Total barcoded reads |
| QC06 | GEX | Reads in passing cells | >= 70.00% | 原生 Fraction of reads in passing cells |
| QC07 | GEX | Median genes/passing cell | >= 1,000 | Passing cells 的基因数中位数 |
| QC08 | GEX | Median molecules/passing cell | >= 2,000 | Passing cells 的分子数中位数 |
| QC09 | CRISPR capture | Barcoded feature reads 比例 | >= 70.00% | Total barcoded feature reads / Total input feature reads |
| QC10 | CRISPR capture | Feature matching rate | >= 60.00% | Feature matching reads / Total barcoded feature reads |
| QC11 | CRISPR capture | Feature reads in passing cells | >= 70.00% | 原生 Fraction of feature reads in passing cells |
| QC12 | CRISPR capture | Passing cells with features 比例 | >= 70.00% | 原生 Fraction passing cells with CRISPR reads；分母为 passing cells |
| QC13 | Guide assignment | Singlet / passing cells | >= 25.00% | 恰好一条阳性 guide 的 passing cells / passing cells |
| QC14 | Guide assignment | Multi-guide / passing cells | <= 20.00% | 至少两条阳性 guides 的 passing cells / passing cells |
| QC15 | Guide assignment | Unassigned / passing cells | <= 60.00% | 无阳性 guide 的 passing cells / passing cells |
| QC16 | NTC | NTC singlet cells | >= 100 | 仅分配一条 NTC guide 的 cells；不能用 feature-positive cells 代替 |
| QC17 | NTC | Active NTC guides | >= 2 | 每条具有至少 5 个 singlet cells 的 NTC guides 数 |
| QC18 | Target coverage | Observed / designed non-NTC targets | >= 90.00% | 有 singlet 的非 NTC targets / 设计中的非 NTC targets |
| QC19 | Target coverage | Eligible / observed targets | >= 25.00% | 至少 20 singlets 的 targets / observed targets |
| QC20 | Guide support | Eligible targets with >=2 active guides | >= 70.00% | 至少两条 active guides 的 eligible targets / eligible targets |
| QC21 | Guide balance | Guide-balanced eligible targets | >= 50.00% | 至少两条 active guides 且 largest active-guide share <=0.5 的 eligible targets / eligible targets |

## 统计单位与取数

- GEX 与 CRISPR input reads 分开统计，R1/R2 不重复计数；以 DRAGEN 对应 input reads 为准。GEX 总 reads 要求为 passing cells × 10,000；CRISPR 总 reads 要求为 passing cells × 800。不得相加后代替单库深度。
- Q30 为 GEX+CRISPR 总体碱基比例，不替换为 GEX-only Q30。Guide assignment 以原生 GMM 阳性分配为准，feature-positive 不等于阳性 guide assignment；guide singlet 不等于已证明不存在生物学 doublet。
- Active guide 为 >=5 singlet cells；eligible target 为 >=20 singlet cells。NTC 由 guide reference 中的 Non-Targeting Control 定义；observed/designed target 均排除 NTC。Largest-guide share = 最大 active guide 的 singlet 数 / 该 target 的全部 singlet 数，沿用现有代码定义，不改为 active-guide-only pool。
- QC21 必须同时满足 eligible、>=2 active guides 和最大 active-guide share <=0.5；不能把仅有一条 active guide 的 target 算作 balanced。
- 原始计数可得时，保存分子、分母及未截断的计算值；原生仅给出比例时保留其实际输出精度，分子记 NA，不从已舍入比例反推整数。
- 共用指标的客户可读显示与现有 DRAGEN 报告精度对齐；精确计数重算值另存 TSV/JSON。应注明因输出精度不同产生的末位差异。
- 缺少投入细胞数时，不推算捕获效率或 reads/input cell。线粒体 reads 比例不能替代每细胞线粒体表达比例。

## 技术与结果一致性

1. 原生 DRAGEN CRISPR 全量运行完成，使用当前样本的运行完成记录与输出；不能用 pilot 或历史结果代替。
2. 自有报告、DRAGEN 官方报告与原生源文件的 35 项共用指标全部匹配；计数严格相等，比例按官方现有输出精度核对。
3. Singlet + multi-guide + unassigned = passing cells，且 barcode 不重复；NTC、guide、target 计数从当前 assignment 重新核对。
4. Evidence 分类数之和 = observed targets；eligible 数、各层 funnel、guide coverage 和 target coverage 与源分配一致；分类规则保持不变。
5. H5AD barcode 顺序与 filtered MEX、assignment 一致，gene feature 顺序一致；已有 X_umap 存在且坐标数等于 passing cells，不重新计算 embedding。
6. 当前报告交付验证全部通过；完整表保留、下载链接有效、核心图齐全、报告离线可用、现有交付文件清单匹配。
7. 报告保留单一样本、有限 NTC、cell-level inference、UNKNOWN perturbation mode、探索性 Reactome、描述性 UMAP 等限制，不能出现超出证据的结论。

## 解释项与适用边界

Sequencing saturation、guide correlation、LOOGO、DE 显著基因数、SUPPORTED 数、Reactome 和 Leiden-cluster sensitivity 须如实报告，但不要求必须出现显著性或 SUPPORTED target 才能通过 QC。QC 通过不等于筛选 hit、编辑成功或机制验证。现有 evidence classification 的阈值和优先级不变。

当前为研究型单样本分析：cell-level p/FDR 不能代替 biological replicate inference；perturbation mode UNKNOWN 时，target RNA 下降或不变均不能证明编辑成功或失败。Reactome 为 exploratory，UMAP 复用既有 X_umap、仅作描述，不能直接推断 perturbation mechanism。

## 未通过时的处理

- 深度未达标：列出按当前 passing cell 数计算的 reads 缺口；补测建议不是补测后必然通过的保证。
- Guide assignment 或覆盖未达标：记录未分配比例、guide 覆盖与 target 覆盖，建议复查富集、文库及实验设计；不得自动归因于某一种实验失败。
- 技术一致性未通过：先修复对应数据或展示问题，再按同一标准复核；不修改原始结果追求通过。
- 缺少必审信息：列明缺失的文件或字段；不得以 NA 充当 PASS。

## 执行与文件

读取当前 native metrics、positive/singlet assignment、guide reference、coverage、H5AD/MEX 和现有报告，形成审核 Markdown、逐项 TSV 与 JSON；记录标准版本、审核日期、源文件、分子/分母、实际值、显示精度、阈值、判定及处理建议。审核文件放入交付的 04_Methods_Reproducibility/review/，脚本不得放入交付目录。

自动审核由正式流程模块执行。每次运行使用配置指定的标准版本；QC 未达标仍保留分析结果，但不标记为验收合格。
