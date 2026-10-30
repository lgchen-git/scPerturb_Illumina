# scPerturb Illumina CRISPR workflow

原生 DRAGEN GEX+CRISPR、下游 perturbation 分析、官方报告、自有离线报告及自动交付审核。项目使用固定版本的本地 Docker，不依赖 GPT、LLM 或 AI API。

## 快速入口

```bash
cd scPerturb_Illumina
pipeline/docker/run_snakemake.sh --help
```

新数据按 [运行手册](pipeline/docs/RUNBOOK.md) 准备独立批次配置后运行：

```bash
pipeline/docker/run_snakemake.sh --config config/PROJECT_YYYYMMDD/run.json --through-dragen --cores 8 --dry-run
pipeline/docker/run_snakemake.sh --config config/PROJECT_YYYYMMDD/run.json --through-dragen --cores 8
# 在 ICM 完成 Perturb-seq 并放回约定导出文件后
pipeline/docker/run_snakemake.sh --config config/PROJECT_YYYYMMDD/run.json --cores 8
```

已完成 DRAGEN 的批次可用 `--prepare-icm` 仅生成本地交接包；不会重跑 DRAGEN 或自动上传到 ICM。详见 [ICM 三级分析运行步骤](pipeline/docs/RUNBOOK.md#icm-perturb-seq-三级分析)。

```bash
pipeline/docker/run_snakemake.sh --config config/run.yaml --prepare-icm
# 在 ICM 完成分析并放回 returned/ 导出文件后
pipeline/docker/run_snakemake.sh --config config/run.yaml --report-only
```

已完成批次只更新报告或审核：

```bash
pipeline/docker/run_snakemake.sh --config config/run.yaml --report-only
pipeline/docker/run_snakemake.sh --config config/run.yaml --review-only
```

本地 `config/reference_path.txt` 和 `config/credential_path.txt` 保存已有参考库及凭证的路径；环境变量 `SC_PERTURB_HG38_ROOT`、`DRAGEN_CREDENTIAL` 可覆盖。凭证内容不复制到项目或交付目录，实际凭证文件需权限 600。授权失败即停止，不使用 fallback 或 pilot 代替全量结果。

首次部署、输入格式、配置示例、断点恢复、运行模式、交付位置及故障处理见 **[RUNBOOK.md](pipeline/docs/RUNBOOK.md)**。已有配置不要重新复制模板覆盖；`run.yaml` 内容实际为 JSON。

## 结构与处理顺序

保留七个目录：`config/`、`inputs/`、`reference/`、`outputs/`、`log/`、`pipeline/`、`tmp/`。新批次通过独立 `results_root/log_root/tmp_root` 隔离；不移动历史结果。

`pipeline/Snakefile` 调度五阶段：

1. `01_input_preflight`：FASTQ、样本和 guide 配对检查。
2. `02_dragen_scrna_crispr`：真实 100,000-read 授权预检及全量原生分析。
3. `03_cell_guide_qc`：细胞、guide 分配和表达空间。
4. `04_primary_reanalysis`：DE、multi-guide consistency、LOOGO 和 Reactome。
5. `05_report`：官方报告、35 项指标对齐、自有 HTML、技术验证及自动审核。

可复用代码在 `pipeline/scperturb/`：`stages/` 编排，`dragen/` 原生运行与官方报告，`analysis/` 统计分析，`reporting/` 绘图/HTML/Sankey/NPG 样式，`validation/` 技术验证与定量审核。正式流程不依赖临时脚本或旧 scRNA 项目路径。

## 审核和交付

审核使用版本化配置 `pipeline/assets/delivery_review_v1.1.json`，规则见 [交付审核标准](pipeline/docs/SC_PERTURB_DELIVERY_REVIEW_CRITERIA.md)。GEX 门槛为 10,000、CRISPR 为 800 reads/passing cell。程序和标准文档使用同一套门槛，分析定义与 evidence classification 不变。

全量分析与报告刷新自动生成审核；`--review-only` 只审核既有结果。报告刷新和独立审核均断网，不读取 license，不重跑上游或 UMAP。技术执行失败保留原交付；QC NOT_PASS/PENDING 保留结果并明确标记，不能作为验收合格结果。

交付位于 `<results_root>/<sample_id>/05_report/04_tool_output/<delivery_name>/`。默认名称 `delivery`，可配置 `OUT_日期_{sample_id}`：

- `01_Report/`：样本名前缀自有 HTML。
- `02_Figures/`：图与 logo。
- `03_Tables/`：完整表及官方指标快照。
- `04_Methods_Reproducibility/`：方法、技术验证及 `review/` 下审核 MD/TSV/JSON。
- `05_DRAGEN_Report/`：官方 HTML。

计数、百分比分母、显示精度及 barcode rank 与官方报告核对。完整结果保留 TSV/TSV.GZ，HTML 离线可用，客户交付目录不放脚本。标准版本变化时保留历史审核。

**不要直接清空 tmp**：部分 primary reanalysis 完成记录仍引用其下的结果/缓存，清理前须核对引用。详见运行手册。

## 环境与验证

当前镜像 `scperturb-crispr:dragen4.5.4-1-v1.0.0`，DRAGEN 4.5.4-1、DRAGEN Reports 4.5.4。优先复用已有镜像；需要重建时运行 `pipeline/docker/build_image.sh`，要求本地基础镜像已存在。修改只读挂载的 Python 流程代码无需重建镜像。

`pipeline/tests/` 覆盖官方指标范围与精度、报告失败时保留旧结果、审核边界和不同样本复用。完整测试在生产容器中运行，不调用 DRAGEN 授权或重跑全量数据。

单个 biological sample；cell-level p/FDR 不代替 biological replicate inference。NTC 支持有限；perturbation mode UNKNOWN；Reactome exploratory；报告复用既有 X_umap。QC 通过不证明编辑成功、筛选 hit 或 perturbation mechanism。
