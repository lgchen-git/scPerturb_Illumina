"""Scoped DRAGEN metric extraction, independent parity, and report-only QC assets."""
from __future__ import annotations
import csv
import html
import json
import re
from pathlib import Path
from scperturb.dragen.reports import parse_number

# group, official table id, official label, source metric, column, percent multiplier, denominator
SPECS = [
 ('gex','scrna_metrics_table_grid','Input Reads','Total input reads',3,1,'GEX input reads'),
 ('gex','scrna_metrics_table_grid','% Barcoded Reads','Total barcoded reads',4,100,'GEX input reads; official rounded fraction'),
 ('gex','scrna_metrics_table_grid','% Mapped to Transcriptome','Mapped reads',4,100,'Barcoded GEX reads'),
 ('gex','scrna_metrics_table_grid','% Sequencing Saturation','Sequencing saturation',3,100,'DRAGEN molecule-count definition; descriptive'),
 ('gex','scrna_metrics_table_grid','Passing Cells','Passing cells',3,1,'DRAGEN passing cells'),
 ('gex','scrna_metrics_table_grid','% Reads in Passing Cells','Fraction of reads in passing cells',3,100,'DRAGEN GEX read denominator'),
 ('gex','scrna_metrics_table_grid','Mean Reads per Cell','Mean reads per cell (Total input reads / Passing cells)',3,1,'GEX input reads / passing cells'),
 ('gex','scrna_metrics_table_grid','Median Molecules per Cell','Median molecules per passing cell',3,1,'Passing cells'),
 ('gex','scrna_metrics_table_grid','Median Genes per Cell','Median genes per passing cell',3,1,'Passing cells'),
 ('gex','scrna_metrics_table_grid','Genes Detected','Total genes detected',3,1,'Genes detected in passing cells; not reference gene features'),
 ('mapping','scrna_mapping_metrics_table_grid','% Mapped to Transcriptome','Mapped reads',4,100,'Barcoded GEX reads'),
 ('mapping','scrna_mapping_metrics_table_grid','% Mapped to Genome','__genome',3,1,'Barcoded GEX reads; excludes non-matching reads'),
 ('mapping','scrna_mapping_metrics_table_grid','% Exon Matching','Unique exon matching reads',4,100,'Barcoded GEX reads'),
 ('mapping','scrna_mapping_metrics_table_grid','% Intron Matching','Unique intron matching reads',4,100,'Barcoded GEX reads'),
 ('mapping','scrna_mapping_metrics_table_grid','% Mitochondrial','Mitochondrial reads',4,100,'Barcoded GEX reads'),
 ('mapping','scrna_mapping_metrics_table_grid','% Antisense','Filtered antisense reads',4,100,'Barcoded GEX reads'),
 ('mapping','scrna_mapping_metrics_table_grid','% Ambiguous','Filtered ambiguously matching reads',4,100,'Barcoded GEX reads'),
 ('mapping','scrna_mapping_metrics_table_grid','% Non-Matching','Filtered non-matching reads',4,100,'Barcoded GEX reads'),
 ('feature','feature_metrics_table_grid','Input Feature Reads','Total input feature reads',3,1,'CRISPR feature input reads'),
 ('feature','feature_metrics_table_grid','% Barcoded Feature Reads','Total barcoded feature reads',4,100,'Feature input reads; official rounded fraction'),
 ('feature','feature_metrics_table_grid','% Feature Matching Reads','Feature matching reads',4,100,'Barcoded feature reads'),
 ('feature','feature_metrics_table_grid','% Feature Reads in Cells','Fraction of feature reads in passing cells',3,100,'DRAGEN feature-read denominator'),
 ('feature','feature_metrics_table_grid','% Feature Matching Reads in Cells','CRISPR fraction valid reads in cells',3,100,'DRAGEN valid feature-read denominator'),
 ('feature','feature_metrics_table_grid','% Cells with Features','Fraction passing cells with CRISPR reads',3,100,'All passing cells'),
 ('feature','feature_metrics_table_grid','Mean Feature Reads per Cell','Mean feature reads per cell (Total input feature reads / Passing cells)',3,1,'Feature input reads / passing cells'),
 ('feature','feature_metrics_table_grid','Median Feature Molecules per Cell','Median feature molecules per passing cell',3,1,'All passing cells; not singlet-only feature UMI median'),
 ('feature','feature_metrics_table_grid','Median Features per Cell','Median features per passing cell',3,1,'All passing cells'),
 ('feature','feature_metrics_table_grid','Features Detected','Total features detected',3,1,'Detected features; not total guide definitions'),
 ('guide','guide_assignment_table_grid','Passing Cells with Fetures','__feature_cells',3,1,'Official guide_metrics.csv: Total cells (guide denominator)'),
 ('guide','guide_assignment_table_grid','Cells with >= 1 Guide','__positive',3,1,'Assigned passing cells'),
 ('guide','guide_assignment_table_grid','% Cells with >= 1 Guide','__positive_pct',3,1,'Official guide_metrics.csv: Total cells'),
 ('guide','guide_assignment_table_grid','Cells with >= 2 Guides','__multi',3,1,'Multi-guide passing cells'),
 ('guide','guide_assignment_table_grid','% Cells with >= 2 Guide','__multi_pct',3,1,'Official guide_metrics.csv: Total cells'),
 ('quality','mapping_qc_table_grid','% Q30','__q30',3,1,'Combined GEX + CRISPR input bases'),
 ('quality','read_metrics_table_grid','Input Reads','__combined_reads',3,1,'Combined GEX + CRISPR reads'),
]


def extract_scoped(report: Path):
    source = report.read_text(encoding='utf-8')
    decoder = json.JSONDecoder()
    out = {}
    for match in re.finditer(r'rowData\s*:\s*\[', source):
        try:
            rows, _ = decoder.raw_decode(source[source.index('[', match.start()):])
        except ValueError:
            continue
        scope = re.findall(r'id="([^"\n]+)"', source[:match.start()])[-1]
        for row in rows:
            if not isinstance(row, dict) or 'metric' not in row or 'value' not in row:
                continue
            key = (scope, row['metric'])
            if key in out:
                if parse_number(str(out[key]['value'])) != parse_number(str(row['value'])):
                    raise ValueError(f'Conflicting duplicate official scoped metric: {key}')
                continue
            out[key] = row
    return out


def extract_rank(report: Path):
    source = report.read_text(encoding='utf-8')
    start = source.index('id="barcode_rank-subsection"')
    start = source.index('var rawTraces = ', start) + len('var rawTraces = ')
    traces, _ = json.JSONDecoder().raw_decode(source[start:])
    if [t['label'] for t in traces] != ['Cells', 'Background']:
        raise ValueError('Unexpected official barcode rank series')
    return traces


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def expected_values(native: Path, sample: str, qc: dict):
    with (native / f'{sample}.scRNA_metrics.csv').open() as handle:
        rows = list(csv.reader(handle))
    values = {row[2]: row for row in rows}
    with (native / f'{sample}.scRNA.guide_metrics.csv').open() as handle:
        guide = {row[2]: int(row[3]) for row in csv.reader(handle)}
    positive = guide['Total cells with one or more guides detected']
    multi = guide['Total cells with two or more guides detected']
    if positive != qc['positive_assignment_cells'] or multi != qc['assignment_status']['multi']:
        raise ValueError('Native guide metrics disagree with downstream assignment counts')
    passing = 0; feature_cells = 0
    with (native / f'{sample}.scRNA.barcodeSummary.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['Filter'] == 'PASS':
                passing += 1; feature_cells += int(int(row['Features']) > 0)
    if passing != qc['filtered_cells'] or not feature_cells:
        raise ValueError(f'Native barcode passing/feature counts invalid: {passing}, {feature_cells}')
    feature_cells = guide['Total cells']  # Official guide denominator; not barcodeSummary Features>0.
    with (native / f'{sample}.mapping_metrics.csv').open() as handle:
        mapping = list(csv.reader(handle))
    q30 = next(float(row[4]) for row in mapping if row[:3] == ['MAPPING/ALIGNING SUMMARY','','Q30 bases'])
    combined = next(int(row[3]) for row in mapping if row[:3] == ['MAPPING/ALIGNING SUMMARY','','Total input reads'])
    n = float(values['Total barcoded reads'][3]); nonmatching = float(values['Filtered non-matching reads'][3])
    extra = {'__genome': (1 - nonmatching/n)*100, '__feature_cells':feature_cells,
             '__positive':positive, '__multi':multi, '__positive_pct':positive/feature_cells*100,
             '__multi_pct':multi/feature_cells*100, '__q30':q30, '__combined_reads':combined}
    return values, extra


def display(value: str) -> str:
    x = parse_number(value)
    if '%' in value: return f'{x:.2f}%'
    if x.is_integer(): return f'{int(x):,}'
    return value


def matches_official(raw, text):
    observed = parse_number(text)
    if '%' not in text and '.' not in text:
        return float(raw).is_integer() and float(raw) == observed
    precision = len(text.replace('%','').split('.')[-1]) if '.' in text else 0
    return round(float(raw), precision) == observed


def prepare_alignment(root: Path, native: Path, qc_path: Path):
    project = json.loads((root/'03_Tables/project_summary.json').read_text())
    sample = project['sample_id']
    report = root/'05_DRAGEN_Report'/f'{sample}_DRAGEN_report.html'
    official = extract_scoped(report)
    qc = json.loads(qc_path.read_text())
    values, extra = expected_values(native, sample, qc)
    snapshot=[]; parity=[]
    for group, scope, label, metric, col, factor, denominator in SPECS:
        row = official[(scope,label)]
        raw = extra[metric] if metric.startswith('__') else float(values[metric][col])*factor
        observed = parse_number(row['value'])
        precision = len(row['value'].replace('%','').split('.')[-1]) if '.' in row['value'] else 0
        if not matches_official(raw,row['value']):
            raise ValueError(f'Official/raw mismatch {scope}/{label}: {observed} vs {raw}')
        percent = '%' in row['value']
        if percent and not 0 <= observed <= 100: raise ValueError(f'Invalid fraction: {label}')
        # Shared project-summary values must also agree with the independent native metrics.
        project_keys = {'Passing Cells':'filtered_cells','Median Molecules per Cell':'median_molecules_per_cell',
                        'Median Genes per Cell':'median_genes_per_cell'}
        if label in project_keys and float(project[project_keys[label]]) != raw:
            raise ValueError(f'Project summary disagrees with native metric: {label}')
        source_suffix = '.scRNA_metrics.csv'
        source_metric = metric
        if metric in {'__positive','__multi','__positive_pct','__multi_pct','__feature_cells'}:
            source_suffix = '.scRNA.guide_metrics.csv'
            source_metric = {'__feature_cells':'Total cells','__positive':'Total cells with one or more guides detected','__multi':'Total cells with two or more guides detected','__positive_pct':'Total cells with one or more guides detected / Total cells','__multi_pct':'Total cells with two or more guides detected / Total cells'}[metric]
        elif metric in {'__q30','__combined_reads'}:
            source_suffix = '.mapping_metrics.csv'
            source_metric = 'MAPPING/ALIGNING SUMMARY: ' + ('Q30 bases, percentage' if metric == '__q30' else 'Total input reads')
        elif metric == '__genome':
            source_metric = '(Total barcoded reads - Filtered non-matching reads) / Total barcoded reads * 100'
        snapshot.append(dict(group=group,scope=scope,metric=label,source_file=sample+source_suffix,source_metric=source_metric,official_value=row['value'],
                             display_value=display(row['value']),source_value=raw,
                             precision=precision,unit='percent' if percent else ('mean/cell' if label.startswith('Mean ') else 'count'),
                             denominator=denominator,description=row.get('desc',row.get('tooltip',''))))
        parity.append(dict(sample_id=sample,metric=label,official_value=row['value'],workflow_value=raw,
                           unit='percent' if percent else 'count',comparison=f'official precision: {precision}',
                           source=f'{sample}{source_suffix}: {source_metric}; scope={scope}',status='PASS'))
    if int(qc['assignment_status']['singlet'])+int(qc['assignment_status']['multi'])+int(qc['assignment_status']['none']) != int(project['filtered_cells']):
        raise ValueError('Invalid passing-cell assignment partition')
    method = root/'04_Methods_Reproducibility/methods/dragen_official_alignment.md'
    method.write_text(
        '# DRAGEN official metric alignment\n\n'
        'Shared metrics are keyed by official table scope and metric name. The two Input Reads rows refer to GEX and combined GEX + CRISPR respectively.\n\n'
        'Counts must match exactly. Percentages are taken from the official reported field, converted once, and displayed to two decimals without inventing precision. Means retain the official precision. Native CSV values are checked independently at official display precision.\n\n'
        'Official guide assignment percentages use Total cells in native .scRNA.guide_metrics.csv as the denominator. The % Cells with Features field comes separately from Fraction passing cells with CRISPR reads in .scRNA_metrics.csv. It must not be used to reconstruct the guide denominator; barcodeSummary Features>0 is also a separate detection statistic. Screen singlet/multi/unassigned fractions use all DRAGEN passing cells.\n\n'
        'Genome mapping, transcriptome matching, MAPQ and base Q30 are different metrics. Main Q30 is the combined-library official percentage. Raw read-flow count ratios remain exact-count calculations and can differ in rounding from percentages supplied by DRAGEN.\n\n'
        'Barcode rank uses the exact official Cells / Background coordinates. Cell calls, threshold, DE, guide effects, LOOGO, evidence classification, Reactome and existing X_umap were not recalculated during this refresh.\n\n'
        'Metric source files and formulas are recorded in 03_Tables/dragen_official_metrics.tsv. Report-only refresh command: pipeline/docker/run_snakemake.sh --config config/run.yaml --report-only.\n')
    write_tsv(root/'03_Tables/dragen_official_metrics.tsv',snapshot)
    write_tsv(root/'04_Methods_Reproducibility/validation/dragen_report_parity.tsv',parity)
    traces=extract_rank(report)
    if traces[0]['data'][-1]['x'] != project['filtered_cells'] or traces[1]['data'][0]['x'] != project['filtered_cells']+1:
        raise ValueError('Barcode rank boundary disagrees with passing cells')
    (root/'03_Tables/dragen_barcode_rank.json').write_text(json.dumps(traces,separators=(',',':'))+'\n')
    render_assets(root,snapshot,traces)
    manifest_path = root/'03_Tables/figure_manifest.tsv'
    with manifest_path.open() as handle: manifest = list(csv.DictReader(handle,delimiter='\t'))
    for row in manifest:
        if row['file'] == 'fig02_gex_crispr_qc.png':
            row['caption'] = 'Official DRAGEN values: reads in passing cells, saturation, genome mapping, combined GEX + CRISPR Q30 and passing-cell medians; scopes and denominators are listed below.'
    manifest = [row for row in manifest if row['file'] != 'fig14_dragen_barcode_rank.png']
    manifest.append(dict(file='fig14_dragen_barcode_rank.png',title='DRAGEN barcode rank',section='qc',caption='Official Cells and Background coordinates; cell calls and molecule threshold unchanged.',source_analysis_version='DRAGEN 4.5.4-1 official report',report_revision='v3.2.1 official-aligned'))
    write_tsv(manifest_path,manifest)
    return snapshot


def render_assets(root, rows, traces):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scperturb.reporting.style import apply_matplotlib_style
    apply_matplotlib_style()
    vals={(r['scope'],r['metric']):parse_number(r['official_value']) for r in rows}
    fig, axes=plt.subplots(1,3,figsize=(12,4.2),layout='constrained')
    panels=[(['GEX reads','Feature reads'],[vals['scrna_metrics_table_grid','% Reads in Passing Cells'],vals['feature_metrics_table_grid','% Feature Reads in Cells']],'Reads in passing cells (%)'),
            (['Saturation','Genome mapped','Q30, combined'],[vals['scrna_metrics_table_grid','% Sequencing Saturation'],vals['scrna_mapping_metrics_table_grid','% Mapped to Genome'],vals['mapping_qc_table_grid','% Q30']],'Official QC (%)'),
            (['Genes/cell','Molecules/cell'],[vals['scrna_metrics_table_grid','Median Genes per Cell'],vals['scrna_metrics_table_grid','Median Molecules per Cell']],'Median per passing cell')]
    for idx,(ax,(labels,numbers,ylabel)) in enumerate(zip(axes,panels)):
        ax.bar(labels,numbers,color=['#4DBBD5','#00A087','#E64B35'][:len(labels)],width=.58)
        ax.set_ylabel(ylabel); ax.spines[['top','right']].set_visible(False)
        if idx<2: ax.set_ylim(0,110)
        else: ax.set_ylim(0,max(numbers)*1.17)
        for i,v in enumerate(numbers): ax.text(i,v+(2 if idx<2 else max(numbers)*.025),f'{v:.2f}%' if idx<2 else f'{int(v):,}',ha='center',fontsize=9)
    fig.savefig(root/'02_Figures/figures/fig02_gex_crispr_qc.png',dpi=220,bbox_inches='tight'); plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.5,4.5),layout='constrained')
    for trace,color in zip(traces,['#00A087','#B8BDC6']):
        ax.plot([p['x'] for p in trace['data']],[p['y'] for p in trace['data']],label=trace['label'],color=color,lw=1.4)
    boundary=traces[0]['data'][-1]
    ax.axvline(boundary['x'],color='#3C5488',lw=.8,ls='--')
    ax.set(xscale='log',yscale='log',xlabel='Barcode rank',ylabel='GEX molecules per barcode')
    ax.set_title(f"Official barcode rank: {boundary['x']:,} cells; molecule threshold {boundary['y']:,}")
    ax.spines[['top','right']].set_visible(False); ax.legend(frameon=False)
    fig.savefig(root/'02_Figures/figures/fig14_dragen_barcode_rank.png',dpi=240,bbox_inches='tight'); plt.close(fig)


def alignment_html(root: Path):
    with (root/'03_Tables/dragen_official_metrics.tsv').open() as handle: rows=list(csv.DictReader(handle,delimiter='\t'))
    def table(groups,table_id):
        items=[r for r in rows if r['group'] in groups]
        body=''.join('<tr data-official-scope="'+html.escape(r['scope'])+'" data-official-metric="'+html.escape(r['metric'],quote=True)+'"><td>'+html.escape(r['metric'].replace('Fetures','Features'))+'</td><td data-official-value>'+html.escape(r['display_value'])+'</td><td>'+html.escape(r['denominator'])+'</td></tr>' for r in items)
        return f'<div class="table-scroll"><table class="data-table" id="{table_id}"><thead><tr><th>Official metric</th><th>Value</th><th>Scope / denominator</th></tr></thead><tbody>{body}</tbody></table></div>'
    qc='<div class="two-col"><article class="panel"><h3>GEX · DRAGEN 官方口径</h3>'+table({'gex'},'official-gex')+'</article><article class="panel"><h3>CRISPR capture · DRAGEN 官方口径</h3>'+table({'feature'},'official-feature')+'</article></div>'
    qc+='<details class="appendix"><summary>官方 mapping 与合并文库质量明细</summary>'+table({'mapping','quality'},'official-mapping')+'<p>Combined input reads = GEX + CRISPR。Q30 为合并文库口径；genome mapping 与 transcriptome matching 不是 MAPQ≥40。</p></details>'
    guide='<h3>官方 Guide assignment（feature-positive passing cells 分母）</h3>'+table({'guide'},'official-guide')+'<p>官方 guide 百分比分母来自原生 guide_metrics.csv 的 Total cells；上方筛选比例使用全部 passing cells。% Cells with Features 另取自 scRNA_metrics.csv 的 CRISPR read-detection 指标，不能用它反推 guide 分母；不同口径的百分比不能直接比较。</p>'
    qc+='<p class="file-note">比例统一为两位小数，保留官方原有精度，不从已舍入比例反推计数。<a href="../03_Tables/dragen_official_metrics.tsv">指标与来源映射 TSV</a> · <a href="../05_DRAGEN_Report/'+html.escape(json.loads((root/'03_Tables/project_summary.json').read_text())['sample_id'])+'_DRAGEN_report.html">DRAGEN 官方完整报告</a></p>'
    qc = '<style>table[id^=official-]{min-width:0;width:100%;table-layout:fixed}table[id^=official-] th:nth-child(1){width:35%}table[id^=official-] th:nth-child(2){width:20%}table[id^=official-] th:nth-child(3){width:45%}table[id^=official-] td:nth-child(2){text-align:right;white-space:nowrap}table[id^=official-] td{overflow-wrap:anywhere}</style>' + qc
    return qc,guide


def verify_alignment(root: Path):
    from html.parser import HTMLParser
    project=json.loads((root/'03_Tables/project_summary.json').read_text())
    official=extract_scoped(root/'05_DRAGEN_Report'/f"{project['sample_id']}_DRAGEN_report.html")
    with (root/'03_Tables/dragen_official_metrics.tsv').open() as h: rows=list(csv.DictReader(h,delimiter='\t'))
    class Parser(HTMLParser):
        def __init__(self): super().__init__(); self.key=None; self.capture=False; self.values={}
        def handle_starttag(self,tag,attrs):
            a=dict(attrs)
            if tag=='tr' and 'data-official-scope' in a: self.key=(a['data-official-scope'],a['data-official-metric'])
            if tag=='td' and 'data-official-value' in a: self.capture=True
        def handle_data(self,data):
            if self.capture and self.key: self.values[self.key]=self.values.get(self.key,'')+data
        def handle_endtag(self,tag):
            if tag=='td': self.capture=False
            if tag=='tr': self.key=None
    p=Parser(); p.feed((root/'01_Report'/f"{project['sample_id']}_scPerturb_CRISPR_report.html").read_text())
    expected={(s,l) for _,s,l,*_ in SPECS}
    if set(p.values)!=expected or {(r['scope'],r['metric']) for r in rows}!=expected: raise ValueError('Missing/extra official metrics in final HTML or table')
    for r in rows:
        k=(r['scope'],r['metric']); reference=official[k]['value']
        if p.values[k]!=display(reference) or r['display_value']!=display(reference): raise ValueError(f'HTML/official mismatch: {k}')
        if not matches_official(float(r['source_value']),reference): raise ValueError(f'Raw-source parity mismatch: {k}')
    saved=json.loads((root/'03_Tables/dragen_barcode_rank.json').read_text())
    html_text=(root/'01_Report'/f"{project['sample_id']}_scPerturb_CRISPR_report.html").read_text()
    embedded=re.search(r"<script type='application/json' id='official-barcode-rank-data'>(.*?)</script>",html_text,re.S)
    if not embedded or json.loads(embedded.group(1)) != saved: raise ValueError('HTML barcode-rank coordinates differ from delivered data')
    if saved!=extract_rank(root/'05_DRAGEN_Report'/f"{project['sample_id']}_DRAGEN_report.html"): raise ValueError('Barcode rank coordinates changed')
    if not (root/'02_Figures/figures/fig14_dragen_barcode_rank.png').is_file(): raise ValueError('Missing barcode rank figure')
    return len(rows)
