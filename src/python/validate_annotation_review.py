"""比较两轮注释及绘图数值；不重新计算或修改生物学结论。"""
import csv
import json
import re
import sys
from pathlib import Path
from pypdf import PdfReader

project, old_dir, new_dir, mapping_path, report_path = sys.argv[1:]
project = Path(project)
old_dir, new_dir = project / old_dir, project / new_dir


def tsv(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def keyed(rows):
    result = {(r['subgroup'], r['subgroup_cluster']): r for r in rows}
    assert len(result) == len(rows)
    return result


old = keyed(tsv(old_dir / 'cluster_annotations.tsv'))
new = keyed(tsv(new_dir / 'cluster_annotations.tsv'))
config = keyed(tsv(project / mapping_path))
assert old.keys() == new.keys() == config.keys()
changes, confidence, functions = [], [], []
for key, row in new.items():
    prior = old[key]
    assert all(row[k] == v for k, v in config[key].items()), key
    assert row['previous_annotation'] == prior['annotation'], key
    for field in ('cells', 'samples_present', 'max_sample_fraction', 'biological_lineage'):
        assert row[field] == prior[field], (key, field)
    item = dict(subgroup=key[0], cluster=key[1], before=prior['annotation'],
                after=row['annotation'], cells=int(row['cells']))
    if row['annotation'] != prior['annotation']:
        changes.append(item)
    if row['confidence'] != prior['confidence']:
        confidence.append(dict(item, before_confidence=prior['confidence'], after_confidence=row['confidence']))
    if row['neutrophil_function'] != prior['neutrophil_function']:
        functions.append(dict(item, before_function=prior['neutrophil_function'], after_function=row['neutrophil_function']))
groups = list(dict.fromkeys(k[0] for k in new))
for group in groups:
    current = tsv(new_dir / group / 'DotPlot_annotation_data.tsv')
    previous = {(r['cluster'], r['gene']): r for r in tsv(old_dir / group / 'DotPlot_annotation_data.tsv')}
    assert len({r['cluster'] for r in current}) == sum(k[0] == group for k in new)
    for r in current:
        key = (r['cluster'], r['gene'])
        if key in previous:
            for column in ('percent_expressed', 'average_expression', 'color_value'):
                assert abs(float(r[column]) - float(previous[key][column])) < 1e-8, (group, key, column)
    names = {r['annotation'] for key, r in new.items() if key[0] == group}
    for stem in ('UMAP_annotated', 'UMAP_sample_split', 'UMAP_group_split', 'DotPlot_annotation_markers'):
        path = project / 'results/06_subcluster' / group / (stem + '.pdf')
        text = '\n'.join(page.extract_text() for page in PdfReader(path).pages)
        assert not re.search(r'\bC\d+\b', text), path
        # R的PDF字体把ASCII连字符映射为Unicode minus；仅规范化提取文字用于比对。
        compact = re.sub(r'\s+', '', text).replace('\u2212', '-')
        assert all(re.sub(r'\s+', '', label) in compact for label in names), (path, names)
for name in ('neutrophil_cluster.tsv', 'neutrophil_cluster_sample.tsv', 'neutrophil_functional_DotPlot_data.tsv'):
    assert tsv(old_dir / name) == tsv(new_dir / name), name
manifest = tsv(new_dir / 'subgroup_annotation_summary.tsv')
complete = tsv(new_dir / 'complete_object_validation.tsv')[0]
assert all(r['readback_verified'] == 'TRUE' and r['RNA_reductions_graphs_preserved'] == 'TRUE' for r in manifest)
assert complete['readback_verified'] == complete['RNA_reductions_graphs_preserved'] == 'TRUE'
cells = sum(int(r['cells']) for r in new.values())
assert cells == int(complete['cells']) == sum(int(r['cells']) for r in manifest)
assert len(tsv(new_dir / 'review_cluster_QC.tsv')) == len(new)
report = dict(annotation_rows=len(new), cells=cells, annotation_name_changes=changes,
              confidence_changes=confidence, function_changes=functions,
              biological_lineage_unchanged=True, biological_lineage_and_function_fields_unchanged=not functions,
              all_common_DotPlot_values_unchanged=True, all_current_figure_labels_match=True,
              duplicate_label_rows_preserved=True, neutrophil_program_tables_unchanged=True,
              QS_readback_and_RNA_reductions_graphs_verified=True)
(project / report_path).write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(report, ensure_ascii=False, indent=2))
print('ANNOTATION_REVIEW_VALIDATED')
