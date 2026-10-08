"""核对逐群文献用途、已执行检索、实际marker比例和Seurat点图数值。"""
import csv
import json
import sys
from pathlib import Path

project, evidence, mapping_file, search_dir, report_file, *additional_search_dirs = sys.argv[1:]
project = Path(project)
evidence = project / evidence


def tsv(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


mapping = tsv(project / mapping_file)
searches = [row for folder in [search_dir, *additional_search_dirs]
            for path in (project / folder).rglob('searches.json')
            for row in json.loads(path.read_text(encoding='utf-8'))]
topics = {row['topic'] for row in searches}
assert len(topics) == len(searches)
assert all(row['query'] and row['source_url'] and row['retrieved_utc']
           and 'count' in row['result'] for row in searches)
records = json.loads((evidence / 'selected_reference_records.json').read_text(encoding='utf-8'))
references = {str(row['pmid']): row for row in records}
observed = {(row['subgroup'], row['cluster'], row['gene']): float(row['percent'])
            for row in tsv(evidence / 'marker_expression.tsv')}
assigned = set()
for row in mapping:
    pmids = set(row['PMIDs'].split(';'))
    identity = set(filter(None, row['identity_PMIDs'].split(';')))
    states = set(filter(None, row['state_PMIDs'].split(';')))
    assert identity and pmids == identity | states
    assert pmids <= references.keys()
    assert set(row['literature_topics'].split(';')) <= topics
    for field in ('identity_criteria_CN', 'alternative_diagnoses_CN', 'reference_scope_CN',
                  'literature_match_level', 'observed_marker_percent'):
        assert row[field], (row['subgroup'], row['subgroup_cluster'], field)
    assert '<U+' not in str(row)
    for item in row['observed_marker_percent'].split(';'):
        gene, percent = item.rsplit(':', 1)
        key = (row['subgroup'], row['subgroup_cluster'], gene)
        assert abs(float(percent.rstrip('%')) - observed[key]) <= 0.051, key
    assigned.update(pmids)

# 61步直接从当前QS用Seurat算出的参考，与63步最终绘图数据逐值一致。
for group in dict.fromkeys(row['subgroup'] for row in mapping):
    actual = {(r['id'], r['features.plot']): r
              for r in tsv(evidence / group / 'Seurat_DotPlot_reference.tsv')}
    plotted = tsv(evidence / group / 'DotPlot_annotation_data.tsv')
    assert {(r['cluster'], r['gene']) for r in plotted} == actual.keys()
    for row in plotted:
        ref = actual[(row['cluster'], row['gene'])]
        for field, original in (('percent_expressed', 'pct.exp'),
                                ('average_expression', 'avg.exp'),
                                ('color_value', 'avg.exp.scaled')):
            assert abs(float(row[field]) - float(ref[original])) < 1e-8

report = dict(clusters=len(mapping), subgroups=len({r['subgroup'] for r in mapping}),
              executed_queries=len(searches), assigned_PMIDs=len(assigned),
              assigned_references_with_full_text=sum(bool(references[p].get('full_text_saved')) for p in assigned),
              all_clusters_have_search_topics_and_reference_roles=True,
              actual_marker_percentages_verified=True, all_DotPlot_values_match_Seurat=True,
              annotation_names=sorted({r['annotation'] for r in mapping}))
(project / report_file).write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(report, ensure_ascii=False, indent=2))
print('IDENTITY_EVIDENCE_VALIDATED')
