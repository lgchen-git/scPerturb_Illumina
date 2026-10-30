"""Reusable quantitative review of completed native scPerturb results."""
from __future__ import annotations
import argparse
import csv
import gzip
import json
import math
import shutil
from collections import Counter
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
from scperturb.common import delivery_name, project_path
from scperturb.validation.delivery import validate, write_manifest
from scperturb.dragen.alignment import verify_alignment, SPECS, extract_scoped, matches_official

DEFAULT_CRITERIA = "pipeline/assets/delivery_review_v1.1.json"


def load_criteria(path: Path) -> dict:
    criteria = json.loads(path.read_text(encoding="utf-8"))
    rules = criteria.get("rules", [])
    if not criteria.get("version") or [r["id"] for r in rules] != [f"QC{i:02d}" for i in range(1,22)]:
        raise ValueError("Review criteria must contain version and ordered QC01..QC21")
    for r in rules:
        value = r["threshold"]
        if r["operator"] not in {">=", "<="} or not isinstance(value,(int,float)) or not math.isfinite(value) or value < 0:
            raise ValueError(f"Invalid review rule: {r['id']}")
        if r["unit"] == "fraction" and value > 1:
            raise ValueError(f"Fraction thresholds must be 0..1: {r['id']}")
    return criteria


def evaluate(value, operator, threshold, denominator=None):
    if value is None or not math.isfinite(value) or denominator == 0:
        return "PENDING"
    return "PASS" if (value >= threshold if operator == ">=" else value <= threshold) else "NOT_PASS"


def overall_result(rows):
    states = {r["status"] for r in rows}
    return "NOT_PASS" if "NOT_PASS" in states else "PENDING" if "PENDING" in states else "PASS"


def fmt(value, unit):
    if value is None or not math.isfinite(value): return "NA"
    if unit == "fraction": return f"{100*value:.2f}%"
    if unit in {"count", "checks"}: return f"{int(value):,}"
    return f"{value:,.2f}"


def criteria_markdown(criteria):
    lines = ["# scPerturb 客户交付审核标准", "",
             f"版本：{criteria['version']}；生效日期：{criteria['effective_date']}。", "",
             "## 判定规则", "",
             "全部必审项通过为 PASS；任一未达标为 NOT_PASS；无未达标但必审信息不足为 PENDING。后两者不能标记为验收合格。分析结果继续保留。",
             "使用源计数或原生最高可用精度判定，边界包含在内；缺失、零分母和非有限值不算通过。改变门槛须更新版本并保留历史审核，不改变分析数据或统计定义。", "",
             "## 定量 QC 门槛", "", "| 编号 | 类别 | 指标 | 通过条件 | 统计口径 |", "|---|---|---|---|---|"]
    for r in criteria["rules"]:
        lines.append(f"| {r['id']} | {r['category']} | {r['metric']} | {r['operator']} {fmt(r['threshold'],r['unit'])} | {r['definition']} |")
    thresholds = {r["id"]:r["threshold"] for r in criteria["rules"]}
    notes = CRITERIA_NOTES.replace("10,000", f"{thresholds['QC01']:,.0f}").replace("× 800", f"× {thresholds['QC02']:,.0f}")
    return ("\n".join(lines) + "\n\n## 统计单位与取数" + notes).rstrip() + "\n"


CRITERIA_NOTES = '\n\n- GEX 与 CRISPR input reads 分开统计，R1/R2 不重复计数；以 DRAGEN 对应 input reads 为准。GEX 总 reads 要求为 passing cells × 10,000；CRISPR 总 reads 要求为 passing cells × 800。不得相加后代替单库深度。\n- Q30 为 GEX+CRISPR 总体碱基比例，不替换为 GEX-only Q30。Guide assignment 以原生 GMM 阳性分配为准，feature-positive 不等于阳性 guide assignment；guide singlet 不等于已证明不存在生物学 doublet。\n- Active guide 为 >=5 singlet cells；eligible target 为 >=20 singlet cells。NTC 由 guide reference 中的 Non-Targeting Control 定义；observed/designed target 均排除 NTC。Largest-guide share = 最大 active guide 的 singlet 数 / 该 target 的全部 singlet 数，沿用现有代码定义，不改为 active-guide-only pool。\n- QC21 必须同时满足 eligible、>=2 active guides 和最大 active-guide share <=0.5；不能把仅有一条 active guide 的 target 算作 balanced。\n- 原始计数可得时，保存分子、分母及未截断的计算值；原生仅给出比例时保留其实际输出精度，分子记 NA，不从已舍入比例反推整数。\n- 共用指标的客户可读显示与现有 DRAGEN 报告精度对齐；精确计数重算值另存 TSV/JSON。应注明因输出精度不同产生的末位差异。\n- 缺少投入细胞数时，不推算捕获效率或 reads/input cell。线粒体 reads 比例不能替代每细胞线粒体表达比例。\n\n## 技术与结果一致性\n\n1. 原生 DRAGEN CRISPR 全量运行完成，使用当前样本的运行完成记录与输出；不能用 pilot 或历史结果代替。\n2. 自有报告、DRAGEN 官方报告与原生源文件的 35 项共用指标全部匹配；计数严格相等，比例按官方现有输出精度核对。\n3. Singlet + multi-guide + unassigned = passing cells，且 barcode 不重复；NTC、guide、target 计数从当前 assignment 重新核对。\n4. Evidence 分类数之和 = observed targets；eligible 数、各层 funnel、guide coverage 和 target coverage 与源分配一致；分类规则保持不变。\n5. H5AD barcode 顺序与 filtered MEX、assignment 一致，gene feature 顺序一致；已有 X_umap 存在且坐标数等于 passing cells，不重新计算 embedding。\n6. 当前报告交付验证全部通过；完整表保留、下载链接有效、核心图齐全、报告离线可用、现有交付文件清单匹配。\n7. 报告保留单一样本、有限 NTC、cell-level inference、UNKNOWN perturbation mode、探索性 Reactome、描述性 UMAP 等限制，不能出现超出证据的结论。\n\n## 解释项与适用边界\n\nSequencing saturation、guide correlation、LOOGO、DE 显著基因数、SUPPORTED 数、Reactome 和 Leiden-cluster sensitivity 须如实报告，但不要求必须出现显著性或 SUPPORTED target 才能通过 QC。QC 通过不等于筛选 hit、编辑成功或机制验证。现有 evidence classification 的阈值和优先级不变。\n\n当前为研究型单样本分析：cell-level p/FDR 不能代替 biological replicate inference；perturbation mode UNKNOWN 时，target RNA 下降或不变均不能证明编辑成功或失败。Reactome 为 exploratory，UMAP 复用既有 X_umap、仅作描述，不能直接推断 perturbation mechanism。\n\n## 未通过时的处理\n\n- 深度未达标：列出按当前 passing cell 数计算的 reads 缺口；补测建议不是补测后必然通过的保证。\n- Guide assignment 或覆盖未达标：记录未分配比例、guide 覆盖与 target 覆盖，建议复查富集、文库及实验设计；不得自动归因于某一种实验失败。\n- 技术一致性未通过：先修复对应数据或展示问题，再按同一标准复核；不修改原始结果追求通过。\n- 缺少必审信息：列明缺失的文件或字段；不得以 NA 充当 PASS。\n\n## 执行与文件\n\n读取当前 native metrics、positive/singlet assignment、guide reference、coverage、H5AD/MEX 和现有报告，形成审核 Markdown、逐项 TSV 与 JSON；记录标准版本、审核日期、源文件、分子/分母、实际值、显示精度、阈值、判定及处理建议。审核文件放入交付的 04_Methods_Reproducibility/review/，脚本不得放入交付目录。\n\n自动审核由正式流程模块执行。每次运行使用配置指定的标准版本；QC 未达标仍保留分析结果，但不标记为验收合格。\n\n'


def collect_review(run: dict, sample: str, root: Path) -> dict:
    import h5py
    repo = Path(run["project_root"])
    batch = Path(run["results_root"]) / sample
    native = batch / "02_dragen_scrna_crispr/04_tool_output/dragen"
    tables = root / "03_Tables"
    analysis = json.loads(project_path(run, run["analysis_config"]).read_text())
    control_gene = analysis["guide_reference"]["control_gene"]
    reanalysis = analysis.get("reanalysis", {})
    if (analysis["qc"]["min_target_cells"], reanalysis.get("guide_min_singlet_cells",5), reanalysis.get("max_guide_fraction",.5)) != (20,5,.5):
        raise ValueError("Review v1.1 requires existing analysis definitions: target>=20, active guide>=5, largest share<=0.5")
    criteria = load_criteria(project_path(run, run.get("delivery_review_config", DEFAULT_CRITERIA)))
    rules = [tuple(r[k] for k in ("id","category","metric","operator","threshold","unit","definition")) for r in criteria["rules"]]
    thresholds = {r["id"]:r["threshold"] for r in criteria["rules"]}
    def tsv(p):
        with p.open() as f: return list(csv.DictReader(f,delimiter='\t'))
    def relative(p):
        if p.is_relative_to(root):
            p = batch / "05_report/04_tool_output" / delivery_name(run, sample) / p.relative_to(root)
        return str(p.relative_to(repo))
    metrics_path = native/f'{sample}.scRNA_metrics.csv'
    with metrics_path.open() as f: metrics = {r[2]:r for r in csv.reader(f)}
    def count(key): return int(metrics[key][3])
    def scalar(key): return float(metrics[key][3])
    project=json.loads((tables/'project_summary.json').read_text())
    if project['sample_id'] != sample: raise ValueError('Report summary belongs to another sample')
    ref_path=project_path(run, run['guide_source'])
    with ref_path.open(encoding='utf-8-sig') as f: ref=list(csv.DictReader(f))
    guide_map={r['sgRNA_id']:r['gene'] for r in ref}
    assert len(guide_map)==len(ref)
    ntc_ids={g for g,t in guide_map.items() if t==control_gene}
    designed={t for t in guide_map.values() if t!=control_gene}
    assignment_path=native/f'{sample}.scRNA.positive_cell_guide_assignments.csv'
    with assignment_path.open() as f: native_assignments=list(csv.DictReader(f))
    with gzip.open(native/f'{sample}.scRNA.filtered.barcodes.tsv.gz','rt') as f: passing_barcodes=[x.strip() for x in f]
    native_by_barcode={r['cell_barcode']:r for r in native_assignments}
    assignment_input_valid=len(native_by_barcode)==len(native_assignments) and set(native_by_barcode)<=set(passing_barcodes)
    missing_native_barcodes=len(set(passing_barcodes)-set(native_by_barcode))
    assignments=[native_by_barcode.get(b,dict(cell_barcode=b,num_features='0',feature_call='',num_transcripts='')) for b in passing_barcodes]
    bc=[r['cell_barcode'] for r in assignments]
    partition=Counter(); guide_singlets=Counter(); target_singlets=Counter(); assignment_valid=assignment_input_valid
    for r in assignments:
        n=int(r['num_features']); ids=r['feature_call'].split('|') if r['feature_call'] else []
        assignment_valid &= len(ids)==n and len(set(ids))==len(ids) and all(g in guide_map for g in ids)
        partition['singlet' if n==1 else 'multi' if n>1 else 'unassigned']+=1
        if n==1:
            guide_singlets[ids[0]]+=1
            if ids[0] not in ntc_ids: target_singlets[guide_map[ids[0]]]+=1
    passing=count('Passing cells'); observed=len(target_singlets)
    eligible={t for t,n in target_singlets.items() if n>=20}
    active_by_target={t:{g:n for g,n in guide_singlets.items() if guide_map[g]==t and n>=5} for t in target_singlets}
    multi_active={t for t in eligible if len(active_by_target[t])>=2}
    balanced={t for t in multi_active if max(active_by_target[t].values())/target_singlets[t]<=.5}
    ntc_singlets=sum(guide_singlets[g] for g in ntc_ids)
    active_ntc=sum(guide_singlets[g]>=5 for g in ntc_ids)
    gex=count('Total input reads'); crispr=count('Total input feature reads')
    with (native/f'{sample}.mapping_metrics.csv').open() as f: mapping=list(csv.reader(f))
    q30=float(next(r[4] for r in mapping if r[:3]==['MAPPING/ALIGNING SUMMARY','','Q30 bases']))/100
    q30_row=next(r for r in mapping if r[:3]==['MAPPING/ALIGNING SUMMARY','','Q30 bases'])
    q30_n=int(q30_row[3]); q30_d=None
    snap=tsv(tables/'dragen_official_metrics.tsv')

    values=[
    (gex,passing,relative(metrics_path),None), (crispr,passing,relative(metrics_path),None),
    (q30_n,q30_d,relative(native/f'{sample}.mapping_metrics.csv'),q30),
    (count('Total barcoded reads'),gex,relative(metrics_path),None),
    (count('Mapped reads'),count('Total barcoded reads'),relative(metrics_path),None),
    (None,None,relative(metrics_path),scalar('Fraction of reads in passing cells')),
    (count('Median genes per passing cell'),None,relative(metrics_path),None),
    (count('Median molecules per passing cell'),None,relative(metrics_path),None),
    (count('Total barcoded feature reads'),crispr,relative(metrics_path),None),
    (count('Feature matching reads'),count('Total barcoded feature reads'),relative(metrics_path),None),
    (None,None,relative(metrics_path),scalar('Fraction of feature reads in passing cells')),
    (None,passing,relative(metrics_path),scalar('Fraction passing cells with CRISPR reads')),
    (partition['singlet'],passing,relative(assignment_path),None),
    (partition['multi'],passing,relative(assignment_path),None),
    (partition['unassigned'],passing,relative(assignment_path),None),
    (ntc_singlets,None,relative(assignment_path)+'; '+relative(ref_path),None),
    (active_ntc,None,relative(assignment_path)+'; '+relative(ref_path),None),
    (observed,len(designed),relative(assignment_path)+'; '+relative(ref_path),None),
    (len(eligible),observed,relative(assignment_path)+'; '+relative(ref_path),None),
    (len(multi_active),len(eligible),relative(assignment_path)+'; '+relative(ref_path),None),
    (len(balanced),len(eligible),relative(assignment_path)+'; '+relative(ref_path),None)]
    # Keep the displayed shared metrics identical to the current official report.
    displays={2:('quality','% Q30'),3:('gex','% Barcoded Reads'),4:('gex','% Mapped to Transcriptome'),5:('gex','% Reads in Passing Cells'),8:('feature','% Barcoded Feature Reads'),9:('feature','% Feature Matching Reads'),10:('feature','% Feature Reads in Cells'),11:('feature','% Cells with Features')}
    checks=[]
    for i,(rule,v) in enumerate(zip(rules,values)):
        rid,group,name,op,threshold,unit,definition=rule
        numerator,denominator,source,raw=v
        value=raw if raw is not None else numerator/denominator if denominator else numerator if denominator is None else None
        status=evaluate(value,op,threshold,denominator)
        display=fmt(value,unit)
        if i in displays:
            group_key,label=displays[i]
            candidates=[r for r in snap if r['group']==group_key and r['metric'].lower()==label.lower()]
            if candidates: display=candidates[0]['display_value']
        checks.append(dict(id=rid,category=group,metric=name,value=value,unit=unit,display_value=display,operator=op,threshold=threshold,threshold_display=op+' '+fmt(threshold,unit),numerator=numerator,denominator=denominator,status=status,definition=definition,source_file=source,precision_note='判定使用源计数比值；共用指标显示保持现有报告精度' if raw is None and unit=='fraction' else '使用源文件可用精度'))

    # Recompute consistency from native assignments and current tables.
    technical=[]
    def tech(name,ok,detail,source): technical.append(dict(id='TECH'+str(len(technical)+1).zfill(2),metric=name,status='PASS' if ok else 'NOT_PASS',detail=detail,source_file=source))
    coverage=tsv(tables/'target_coverage_full.tsv'); cov={r['target_gene']:r for r in coverage}
    guides=tsv(tables/'guide_coverage.tsv'); conf=tsv(tables/'target_confidence.tsv')
    tech('Assignment 分区及 barcode 唯一性',len(bc)==passing==sum(partition.values()) and len(set(bc))==passing and assignment_valid,str(dict(partition))+f'；原生表 {len(native_assignments)} 行，按 filtered MEX 全部 passing barcodes 对齐，未列出的 {missing_native_barcodes} 个无阳性 guide',relative(assignment_path))
    tech('Guide coverage 与原生 singlet 一致',len(guides)==len(ref) and {r['sgRNA_id'] for r in guides}==set(guide_map) and all(int(r['singlet_cells'])==guide_singlets[r['sgRNA_id']] and r['target_gene']==guide_map[r['sgRNA_id']] for r in guides),'逐 guide 核对 '+str(len(guides))+' 行',relative(tables/'guide_coverage.tsv'))
    tech('Target coverage 与原生 singlet 一致',set(cov)==set(target_singlets) and all(int(cov[t]['singlet_cells'])==n and int(cov[t]['active_guides_ge5'])==len(active_by_target[t]) for t,n in target_singlets.items()),f'observed={observed}; eligible={len(eligible)}',relative(tables/'target_coverage_full.tsv'))
    funnel_counts=[int(r['count']) for r in tsv(tables/'screen_target_funnel.tsv')]
    consistency_count=sum(t in balanced and float(cov[t]['median_pairwise_spearman'] or 'nan')>=0 and float(cov[t]['median_directional_concordance'] or 'nan')>=.5 for t in cov)
    tech('Target funnel 一致',funnel_counts==[observed,len(eligible),len(multi_active),len(balanced),consistency_count],str(funnel_counts),relative(tables/'screen_target_funnel.tsv'))
    tech('Evidence target 名称及分类分区一致',len(conf)==observed and {r['target_gene'] for r in conf}==set(cov) and sum(project['confidence_counts'].values())==observed and Counter(r.get('confidence_class',r.get('evidence_class','')) for r in conf)==Counter(project['confidence_counts']),str(project['confidence_counts']),relative(tables/'target_confidence.tsv'))
    tech('NTC 计数一致',ntc_singlets==project['ntc_qc']['ntc_singlet_cells'] and active_ntc==project['ntc_qc']['active_ntc_guides_ge5_singlets'] and len(ntc_ids)==project['ntc_qc']['ntc_guide_features'],f'NTC guides={len(ntc_ids)}; active={active_ntc}; singlets={ntc_singlets}',relative(tables/'ntc_guide_qc.tsv'))
    with (native/f'{sample}.scRNA.singlet_cell_guide_assignments.csv').open() as f: native_singlets=list(csv.DictReader(f))
    expected_singlets={r['cell_barcode']:r['feature_call'] for r in assignments if int(r['num_features'])==1}
    tech('原生 singlet 文件与 positive assignment 一致',len(native_singlets)==len(expected_singlets) and {r['cell_barcode']:r['feature_call'] for r in native_singlets}==expected_singlets,'原生两份 assignment 交叉核对',relative(native/f'{sample}.scRNA.singlet_cell_guide_assignments.csv'))
    with gzip.open(native/f'{sample}.scRNA.filtered.barcodes.tsv.gz','rt') as f: mex_bc=[x.strip() for x in f]
    with gzip.open(native/f'{sample}.scRNA.features.tsv.gz','rt') as f: mex_gene=[r[0] for r in csv.reader(f,delimiter='\t') if r[2]=='Gene Expression']
    h5_path=batch/'03_cell_guide_qc/04_tool_output'/f'{sample}_gex_filtered.h5ad'
    with h5py.File(h5_path,'r') as f:
        obs=f['obs'][f['obs'].attrs['_index']].asstr()[:].tolist()
        var=f['var'][f['var'].attrs['_index']].asstr()[:].tolist()
        umap_shape=list(f['obsm/X_umap'].shape)
    tech('H5AD/MEX/assignment barcode 和 feature 顺序',obs==mex_bc==bc and var==mex_gene,f'cells={len(obs)}; genes={len(var)}；本轮直接读取 H5AD metadata 核验',relative(h5_path))
    tech('既有 X_umap 存在',umap_shape==[passing,2],f'X_umap shape={umap_shape}；未计算 embedding',relative(h5_path))
    marker=batch/'02_dragen_scrna_crispr/05_metadata/stage_complete.json'
    completion=json.loads(marker.read_text())
    tech('原生全量 DRAGEN CRISPR 完成',completion.get('status')=='complete' and completion.get('sample_id')==sample and project['dragen_native_crispr'] and project['dragen_finished_normally'] and project['dragen_license_success_in_source_log'] and project['gex_reads']==gex and project['crispr_feature_reads']==crispr,'当前全量完成记录、现有原生完成证据和当前 metrics；审核不调用 license',relative(marker))
    # Independent source -> official parity, without re-running a pipeline stage.
    with (native/f'{sample}.scRNA.guide_metrics.csv').open() as f: gm={r[2]:int(r[3]) for r in csv.reader(f)}
    combined=int(next(r[3] for r in mapping if r[:3]==['MAPPING/ALIGNING SUMMARY','','Total input reads']))
    extra={'__genome':(1-count('Filtered non-matching reads')/count('Total barcoded reads'))*100,'__feature_cells':gm['Total cells'],'__positive':gm['Total cells with one or more guides detected'],'__multi':gm['Total cells with two or more guides detected'],'__positive_pct':gm['Total cells with one or more guides detected']/gm['Total cells']*100,'__multi_pct':gm['Total cells with two or more guides detected']/gm['Total cells']*100,'__q30':q30*100,'__combined_reads':combined}
    official_rows=extract_scoped(root/'05_DRAGEN_Report'/f'{sample}_DRAGEN_report.html')
    parity=[]
    for group,scope,label,key,col,factor,denom in SPECS:
        raw=extra[key] if key.startswith('__') else float(metrics[key][col])*factor
        parity.append(matches_official(raw,official_rows[(scope,label)]['value']))
    tech('35 项原生源文件与官方报告一致',len(parity)==35 and all(parity) and gm['Total cells with one or more guides detected']==partition['singlet']+partition['multi'] and gm['Total cells with two or more guides detected']==partition['multi'] and combined==gex+crispr,f'{sum(parity)}/{len(parity)}；按现有官方显示精度比较',relative(metrics_path))
    tech('35 项自有 HTML/官方 HTML/表格一致',verify_alignment(root)==35,'同时检查 barcode rank coordinates',relative(root/'01_Report'/f'{sample}_scPerturb_CRISPR_report.html'))

    observations={'passing_cells':passing,'gex_reads':gex,'crispr_reads':crispr,'combined_reads':combined,'designed_non_ntc_targets':len(designed),'observed_targets':observed,'eligible_targets':len(eligible),'eligible_with_two_active_guides':len(multi_active),'balanced_eligible_targets':len(balanced),'consistency_gate_targets':consistency_count,'ntc_singlet_cells':ntc_singlets,'active_ntc_guides':active_ntc,'evidence_class_counts':project['confidence_counts'],'sequencing_saturation':project['sequencing_saturation'],'gsea_primary_tie_median':project['gsea_primary_tie_median'],'native_assignment_rows':len(native_assignments),'passing_barcodes_not_listed_in_native_assignments':missing_native_barcodes,'h5ad_barcode_order_verified_this_review':obs==mex_bc==bc,'h5ad_feature_order_verified_this_review':var==mex_gene,'umap_shape':umap_shape,'gex_minimum_input_reads_at_current_passing_cells':passing*thresholds['QC01'],'crispr_minimum_input_reads_at_current_passing_cells':passing*thresholds['QC02'],'gex_read_deficit_at_fixed_cells':max(0,passing*thresholds['QC01']-gex),'crispr_read_deficit_at_fixed_cells':max(0,passing*thresholds['QC02']-crispr)}

    technical.append(dict(id="TECH13", metric="现有报告交付验证", status="PENDING", detail="技术验证尚未执行", source_file="04_Methods_Reproducibility/validation/report_v3_2_1_validation.json"))
    return {"criteria_version":criteria["version"], "criteria":criteria, "sample_id":sample,
            "reviewed_at":datetime.now(ZoneInfo("Asia/Shanghai")).isoformat(timespec="seconds"),
            "checks":checks,"technical_checks":technical,"observations":observations,"analysis_recomputed":False}


def apply_validation(result, validation):
    row = result["technical_checks"][-1]
    passed = sum(r["status"] == "PASS" for r in validation["checks"])
    checks = validation.get("checks", [])
    all_checks_pass = bool(checks) and passed == len(checks)
    row.update(status="PASS" if validation.get("result") == "PASS" and all_checks_pass else "NOT_PASS",
               detail=f"{passed}/{len(checks)} 项技术交付验证通过")
    return result


def write_review(root: Path, result: dict):
    review = root / "04_Methods_Reproducibility/review"
    review.mkdir(parents=True, exist_ok=True)
    base = review / f"{result['sample_id']}_DELIVERY_REVIEW"
    previous = base.with_suffix(".json")
    if previous.exists():
        old = json.loads(previous.read_text())
        if old.get("criteria_version") != result["criteria_version"]:
            history = review / "history" / f"v{old['criteria_version']}"
            history.mkdir(parents=True, exist_ok=True)
            for p in [base.with_suffix(s) for s in (".md",".tsv",".json")] + [review/"SC_PERTURB_DELIVERY_REVIEW_CRITERIA.md"]:
                if p.exists() and not (history/p.name).exists(): shutil.copy2(p,history/p.name)
    all_rows = result["checks"] + result["technical_checks"]
    result.update(overall_result=overall_result(all_rows),
                  quantitative_qc_counts=dict(Counter(r["status"] for r in result["checks"])),
                  technical_counts=dict(Counter(r["status"] for r in result["technical_checks"])),
                  status_counts=dict(Counter(r["status"] for r in all_rows)))
    result["release_allowed"] = result["overall_result"] == "PASS"
    result["criteria_file"] = "SC_PERTURB_DELIVERY_REVIEW_CRITERIA.md"
    (review/result["criteria_file"]).write_text(criteria_markdown(result["criteria"]), encoding="utf-8")
    base.with_suffix(".json").write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    fields = ["id","category","metric","value","unit","display_value","operator","threshold","threshold_display","numerator","denominator","status","definition","source_file","precision_note","detail"]
    with base.with_suffix(".tsv").open("w",encoding="utf-8",newline="") as f:
        writer = csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",restval="NA")
        writer.writeheader()
        for row in all_rows: writer.writerow({k:("NA" if v is None else v) for k,v in row.items()})
    qc = result["checks"]; technical=result["technical_checks"]; o=result["observations"]
    md = [f"# {result['sample_id']} scPerturb 交付审核结果", "",
          f"审核时间：{result['reviewed_at']}；审核标准：v{result['criteria_version']}。", "",
          f"## 总结：{result['overall_result']}", "",
          f"定量 QC：{sum(r['status']=='PASS' for r in qc)}/{len(qc)} 通过；技术与一致性：{sum(r['status']=='PASS' for r in technical)}/{len(technical)} 通过。", "",
          "满足本版验收标准；QC 通过不表示生物学证据增强。" if result["release_allowed"] else "不能标记为验收合格交付；分析结果保留。技术完成与定量 QC 判定分别记录。", "",
          "## 定量 QC", "", "| 编号 | 指标 | 实测值 | 通过条件 | 结果 |", "|---|---|---|---|---|"]
    for r in qc: md.append(f"| {r['id']} | {r['metric']} | {r['display_value']} | {r['threshold_display']} | {r['status']} |")
    md += ["", "共用指标显示与现有 DRAGEN 报告精度一致；TSV/JSON 的 value 保存最高可用精度，判定使用 value，不从已舍入百分比反推计数。", "",
           "## 技术与一致性", "", "| 检查 | 结果 | 说明 |", "|---|---|---|"]
    for r in technical: md.append(f"| {r['metric']} | {r['status']} | {r['detail'].replace('|',' / ')} |")
    md += ["", "## 解读与后续处理", "",
           f"- Passing cells {o['passing_cells']:,}；GEX reads {o['gex_reads']:,}；CRISPR reads {o['crispr_reads']:,}。",
           f"- 按当前细胞数及本版门槛，GEX 需至少 {o['gex_minimum_input_reads_at_current_passing_cells']:,.0f} reads，缺口 {o['gex_read_deficit_at_fixed_cells']:,.0f}；CRISPR 需至少 {o['crispr_minimum_input_reads_at_current_passing_cells']:,.0f}，缺口 {o['crispr_read_deficit_at_fixed_cells']:,.0f}。补测后细胞数可能变化，补足算术缺口不保证改善分配和覆盖。",
           f"- Observed targets {o['observed_targets']:,}，eligible {o['eligible_targets']:,}；SUPPORTED {o['evidence_class_counts'].get('SUPPORTED',0)}，PARTIAL_SUPPORT {o['evidence_class_counts'].get('PARTIAL_SUPPORT',0)}。数量不足或多 guide 不一致时不能提升证据等级。",
           f"- NTC singlets {o['ntc_singlet_cells']:,}，active NTC guides {o['active_ntc_guides']}；达到数量门槛不等于各 cell states 都有充分对照。",
           f"- Sequencing saturation {100*o['sequencing_saturation']:.2f}%；Reactome median ranking tie fraction {100*o['gsea_primary_tie_median']:.2f}%，不作为独立通过线。",
           "", "仅一个 biological sample；cell-level p/FDR 不能替代 biological replicate inference。Perturbation mode UNKNOWN，target RNA 下降或不变不能证明编辑成功或失败。Reactome 为 exploratory；UMAP 使用既有 X_umap、仅作 descriptive 展示，不能直接证明 perturbation mechanism。", "",
           "## 文件", "", f"- [审核标准]({result['criteria_file']})", f"- [逐项 TSV]({base.name}.tsv)", f"- [机器可读 JSON]({base.name}.json)",
           "- 源文件路径以项目根目录为基准；分子、分母、未截断值及精度说明见 TSV/JSON。",
           "- 不重跑分析或 UMAP；脚本不进入交付目录。历史标准变更审核保留在 history/（如有）。"]
    base.with_suffix(".md").write_text("\n".join(md)+"\n",encoding="utf-8")


def prepare_review(run, sample, root):
    result = collect_review(run, sample, root)
    write_review(root,result)
    return result


def finish_review(root, result, validation):
    apply_validation(result,validation)
    write_review(root,result)
    write_manifest(root)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",type=Path,required=True)
    parser.add_argument("--sample-id",required=True)
    args = parser.parse_args()
    run = json.loads(args.config.read_text())
    root = Path(run["results_root"])/args.sample_id/"05_report/04_tool_output"/delivery_name(run,args.sample_id)
    result = prepare_review(run,args.sample_id,root)
    write_manifest(root)
    validation = validate(root,args.sample_id)
    result = finish_review(root,result,validation)
    print(json.dumps({"sample_id":args.sample_id,"delivery_review":result["overall_result"],"quantitative":result["quantitative_qc_counts"],"technical":result["technical_counts"]},ensure_ascii=False))
    return 0 if validation["result"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
