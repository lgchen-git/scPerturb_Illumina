# scPerturb 运行手册

分析在固定版本本地 Docker 容器中执行，不需要 GPT、LLM 或 AI API。原生 DRAGEN 授权仍需厂商服务可用及有效凭证。

## 1. 当前机器直接运行

```bash
cd scPerturb_Illumina
pipeline/docker/run_snakemake.sh --config config/run.yaml --cores 8
```

上述命令使用当前本地配置。**新数据先按第 3 节创建新批次配置**，不要把新 FASTQ 填进已经完成的同一批次输出路径。当前已完成批次只需更新报告或审核时，使用第 4 节命令。

当前环境已准备好：

- 镜像：`scperturb-crispr:dragen4.5.4-1-v1.0.0`；DRAGEN 4.5.4-1、DRAGEN Reports 4.5.4。
- 参考库路径读取本地 `config/reference_path.txt`；该文件只应包含当前机器上的 hg38 参考库路径。
- 凭证路径读取本地 `config/credential_path.txt`，复用既有成功凭证。凭证内容不复制到项目或交付目录；实际凭证文件权限为 `600`。授权过期或额度不足仍会停止原生分析。
- 环境变量 `SC_PERTURB_HG38_ROOT`、`DRAGEN_CREDENTIAL` 可分别覆盖上述路径。
- 最新审核配置：`pipeline/assets/delivery_review_v1.1.json`；GEX 深度 >=10,000、CRISPR 深度 >=800 reads/passing cell，其余要求见[审核标准](SC_PERTURB_DELIVERY_REVIEW_CRITERIA.md)。

## 2. 首次部署或环境恢复

需要 Docker、宿主机 Python 3、Bash，以及本地镜像与参考库。普通运行不需要宿主机安装分析依赖。

```bash
cd scPerturb_Illumina
pipeline/docker/run_snakemake.sh --help
docker image inspect scperturb-crispr:dragen4.5.4-1-v1.0.0
```

优先使用已有镜像。确需重建且本地基础镜像已具备时运行 `pipeline/docker/build_image.sh`；构建复用已有 base-image lock，不自动拉取基础镜像。Python 流程代码通过只读挂载加载，修改流程代码无需重建镜像。

在新机器创建 `config/reference_path.txt`、`config/credential_path.txt`，每个文件只写对应资源的绝对路径。hg38 根目录需要包含运行配置中的 DRAGEN hash、FASTA、GTF 和 STAR reference。Reactome GMT 位于项目 `reference/pathways/`，文件路径由分析配置指定。

本地路径指针、样本配置、凭证、数据和运行产物不提交 Git；可复用代码、配置模板、门槛和文档由 Git 管理。

## 3. 新批次准备

每个样本需四个配对文件：GEX R1/R2、CRISPR R1/R2。Guide CSV 必须具有 `gene,sgRNA_id,sequence` 三列；NTC 默认命名为 `Non-Targeting Control`。每个样本一行，不允许重复 sample_id。当前样本表不接受把同一样本的多个 lane 写成重复行；须先准备好正确配对的输入，不能自动混入旧批次 reads。

目录示例：

```text
inputs/PROJECT_YYYYMMDD/
  guides.csv
  SAMPLE_001_GEX_R1.fastq.gz
  SAMPLE_001_GEX_R2.fastq.gz
  SAMPLE_001_CRISPR_R1.fastq.gz
  SAMPLE_001_CRISPR_R2.fastq.gz
config/PROJECT_YYYYMMDD/
  run.json
  samples.csv
  perturbation.json
```

从模板新建配置，示例不会覆盖当前 `config/run.yaml`：

```bash
cd scPerturb_Illumina
task_batch=PROJECT_YYYYMMDD
mkdir "config/$task_batch"
cp config/run.template.json "config/$task_batch/run.json"
cp config/samples.template.csv "config/$task_batch/samples.csv"
cp config/perturbation.template.json "config/$task_batch/perturbation.json"
python3 - "$task_batch" <<'PY'
import json, sys
from pathlib import Path
batch = sys.argv[1]
p = Path('config') / batch / 'run.json'
c = json.loads(p.read_text())
c.update(project_id=batch,
         sample_sheet=f'/project/config/{batch}/samples.csv',
         analysis_config=f'/project/config/{batch}/perturbation.json',
         guide_source=f'/project/inputs/{batch}/guides.csv',
         results_root=f'/project/outputs/{batch}',
         log_root=f'/project/log/{batch}',
         tmp_root=f'/project/tmp/{batch}',
         delivery_name='OUT_{sample_id}')
p.write_text(json.dumps(c, ensure_ascii=False, indent=2) + '\n')
PY
```

编辑新样本表的 `sample_id,project_id,lane,rglb,gex_r1,gex_r2,crispr_r1,crispr_r2`。FASTQ 路径使用项目相对路径或 `/project/inputs/...`；例如 `inputs/PROJECT_YYYYMMDD/SAMPLE_001_GEX_R1.fastq.gz`。`lane` 使用实际 lane 编号，`rglb` 表示对应样本/文库。`sample_id` 用作文件名和目录名。

运行配置在容器中使用 `/project/...`，参考库挂载在 `/reference/hg38`，不要写宿主机的结果路径。`run.yaml` 历史文件名虽为 YAML，内容仍必须是 JSON。`--config` 参数使用项目相对路径。资源默认 8 线程、32 GiB，按机器条件修改 `resources`；`--cores` 为 Snakemake 总线程预算。

`delivery_name` 可为 `OUT_日期_{sample_id}`。审核规则由 `delivery_review_config` 指向版本化 JSON；缺少该字段的旧配置默认使用 v1.1。更新交付门槛时创建新版本，不要覆盖历史标准。分析的 active guide、eligible target 和 evidence 定义保持现有默认值；当前 v1.1 审核要求 target >=20 singlets、active guide >=5 singlets、largest-guide share <=0.5。

## 4. 运行模式

### 检查新批次 DAG

```bash
pipeline/docker/run_snakemake.sh --config config/PROJECT_YYYYMMDD/run.json --cores 8 --dry-run
```

检查本地输入、参考、凭证文件和 DAG；不运行 DRAGEN、不在线测试授权。日志和必要的运行临时目录仍会建立。

### 全量分析

```bash
pipeline/docker/run_snakemake.sh --config config/PROJECT_YYYYMMDD/run.json --cores 8
```

顺序为：输入检查 → 原生 DRAGEN GEX+CRISPR → cell/guide QC → DE/guide consistency/LOOGO/Reactome → 官方报告、自有报告和自动审核。全量模式已包含真实 100,000-read 授权预检；预检不是全量分析结果，失败时停止，不使用 fallback 代替原生结果。

标准全量命令用于 ICM 已完成并放回返回文件后的本地下游与报告。新样本先按下方 ICM 流程运行到 stage 02 并手动完成 ICM；已有 DRAGEN 输出则先运行 `--prepare-icm`。

### ICM Perturb-seq 三级分析

新样本先运行到 DRAGEN 完成并生成本地交接包：

```bash
pipeline/docker/run_snakemake.sh --config config/PROJECT_YYYYMMDD/run.json --through-dragen --cores 8
```

此模式运行输入预检与 DRAGEN，然后停止在 stage 02；它不会运行本地下游 QC、报告或交付。已经完成 DRAGEN 的样本可只生成交接包，不会重跑 DRAGEN：

```bash
pipeline/docker/run_snakemake.sh --config config/run.yaml --prepare-icm
```

交接包位于 `<results_root>/<sample_id>/02_dragen_scrna_crispr/04_tool_output/icm_perturbseq_handoff/`。将 `input/` 下五个文件手动导入 Illumina Connected Multiomics，选择 **Default: Illumina Single Cell Transcriptomics Perturb-seq**。流程不自动连接或上传到 ICM。完成官方分析后，将表达/细胞组成概览 PNG、guide capture 概览 PNG 和分析摘要 PDF 按 `returned/icm_export_manifest.tsv` 指定的文件名放入 `returned/`。

返回文件到位后再运行全量分析（新样本）或报告刷新（已有完成分析的样本）。报告以 ICM 静态导出替换本地常规 GEX/CRISPR QC 图和 UMAP 探索性图；官方 DRAGEN QC、guide 一致性、LOOGO、靶基因证据分类和相关结构化结果仍保留。缺少 ICM 返回图或摘要时，报告刷新会失败并保留当前正式交付。

### 只诊断授权

```bash
pipeline/docker/run_snakemake.sh --config config/PROJECT_YYYYMMDD/run.json --pilot-only --cores 8
```

结果隔离在 `tmp/native_license_pilot_时间戳/`，不覆盖正式分析结果。该模式会真实调用授权，只在需要诊断时使用，不必在每次全量运行前重复执行。

### 只刷新报告

```bash
pipeline/docker/run_snakemake.sh --config config/run.yaml --report-only --dry-run
pipeline/docker/run_snakemake.sh --config config/run.yaml --report-only
```

从已完成交付和原生输出刷新报告，自动更新审核。必须先放回完整的 ICM PNG/PDF 导出；否则在启动容器前失败。容器断网，不读取 license，不调度上游分析，不重新计算 DE、Reactome、clustering 或 UMAP。技术验证通过后才替换正式交付；失败保留原交付，历史交付保留在该样本 `tmp_root` 下。

### 只重新审核

```bash
pipeline/docker/run_snakemake.sh --config config/run.yaml --review-only --dry-run
pipeline/docker/run_snakemake.sh --config config/run.yaml --review-only
```

只读取现有 native metrics、assignment、guide reference、coverage、H5AD/MEX 和报告，重写审核文件及现有技术验证/文件清单；不改自有 HTML 或分析结果。容器断网，不读取 license。需要已有完整报告、guide CSV、分析配置、native 输出、QC H5AD 和完整 ICM 返回文件。

专用模式 `--through-dragen`、`--prepare-icm`、`--pilot-only`、`--report-only`、`--review-only` 互斥。程序执行成功与 QC 达标分开：QC 为 NOT_PASS/PENDING 时分析结果仍保留，输出明确审核状态，不能据 RUN_COMPLETE 单独认定可验收交付；源文件缺失、解析或技术验证失败按执行错误处理。

### 断点恢复

排除错误后，用**同一个批次配置及同一条全量命令**重新运行，Snakemake 沿用阶段完成记录并恢复未完成任务。不手工伪造完成标记，不删除全部输出来恢复。

完成批次在代码迁移后可能被 Snakemake 判定需重跑；如果仅更新报告或审核，使用对应专用模式，不直接重跑旧批次全量分析。新批次必须使用独立 `results_root/log_root/tmp_root`。

## 5. 如何确认完成及交付位置

五阶段均有 `05_metadata/stage_complete.json`，运行入口执行最终完成检查。交付位置为：

```text
<results_root>/<sample_id>/05_report/04_tool_output/<delivery_name>/
  01_Report/<sample_id>_scPerturb_CRISPR_report.html
  02_Figures/
  03_Tables/
  04_Methods_Reproducibility/
    review/<sample_id>_DELIVERY_REVIEW.md
    review/<sample_id>_DELIVERY_REVIEW.tsv
    review/<sample_id>_DELIVERY_REVIEW.json
    review/SC_PERTURB_DELIVERY_REVIEW_CRITERIA.md
  05_DRAGEN_Report/<sample_id>_DRAGEN_report.html
```

客户交付目录不含脚本。完整分析表保留 TSV/TSV.GZ，报告为离线 HTML；ICM PNG 内嵌于 HTML，官方 ICM 分析摘要 PDF 作为本地文件链接。共用指标对齐 DRAGEN 官方报告的 35 项指标和各自分母；计数精确一致，比例按现有官方输出精度显示。

审核 JSON 的 `overall_result` 是 PASS、NOT_PASS 或 PENDING；只有 PASS 才有 `release_allowed=true`。标准版本变化时保留旧审核到 `review/history/v版本/`。审核 PASS 不等于出现 SUPPORTED target 或证明机制。

## 6. 日志、错误与临时文件

- 总日志：`log/snakemake_时间戳.log`；阶段日志在配置的 `log_root/<sample_id>/`。
- 授权失败：检查有效凭证和原生授权日志。文件存在不代表授权成功，不重复使用失败结果当作全量结果。
- 路径缺失：核对宿主路径指针及容器 `/project`、`/reference/hg38` 路径；不要把两种路径混写。
- 输入或 guide 错误：核对样本库配对、CSV 列、唯一 guide ID、NTC 定义及实际 R1/R2，不通过放宽 QC 修复配对错误。
- 内存不足：先查看失败阶段日志，再调整资源和总 cores；不要改统计方法来绕过错误。
- 报告/审核失败：先核对既有源文件和完成标记，不自动触发上游重算。

**不要直接清空整个 tmp。** `04_primary_reanalysis/05_metadata/reanalysis_root.json` 引用的 `v2_root` 目前位于配置 `tmp_root` 下，阶段完成记录和恢复仍可能依赖它。`report_build` 的历史交付用于回滚；`primary_runtime`、`v3_runtime` 和分析缓存可能用于恢复。需清理时先检查当前完成记录引用，仅清理确认不再需要的独立旧任务目录。不要移动或删除正式交付及原生输出。

## 7. 测试与科学解释

```bash
mkdir -p tmp/tests
TMPDIR="$PWD/tmp/tests" PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=pipeline \
  python3 -m unittest discover -s pipeline/tests -v
```

宿主机缺少分析依赖时，直接在既有镜像中测试：

```bash
mkdir -p tmp/tests
docker run --rm --network none --read-only --user "$(id -u):$(id -g)" \
  --env PYTHONDONTWRITEBYTECODE=1 --env PYTHONPATH=/project/pipeline \
  --env TMPDIR=/project/tmp/tests \
  --entrypoint /opt/conda/envs/illumina_scrna/bin/python \
  --volume "$PWD:/project:rw" --workdir /project \
  scperturb-crispr:dragen4.5.4-1-v1.0.0 \
  -m unittest discover -s pipeline/tests -v
```

完整测试需生产容器的分析依赖（包括 h5py）；宿主机缺依赖时在现有生产容器内执行，临时路径仍指向 `/project/tmp/tests`。测试和 dry-run 不代表新的全量 DRAGEN 分析已完成。

本流程为单样本分析；cell-level p/FDR 不代替 biological replicate inference。NTC 数量和 cell-state 覆盖影响解释。Perturbation mode UNKNOWN，target RNA 变化不能证明编辑成功或失败；Reactome exploratory；报告 UMAP 复用既有 X_umap、仅作 descriptive 展示，不直接证明 perturbation mechanism。
