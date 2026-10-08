"""直接合并各大类正式PDF；原页面保持矢量，排除重绘保留亚群。"""
import csv
import sys
from pathlib import Path
from pypdf import PdfReader, PdfWriter

project = Path(sys.argv[1])
with (project / sys.argv[2]).open(encoding='utf-8-sig', newline='') as fh:
    groups = [row['subgroup'] for row in csv.DictReader(fh, delimiter='\t')]
assert len(groups) == len(set(groups))
output = project / 'results/06_subcluster'
for stem, name in [('UMAP_annotated', 'All_subgroups_UMAP_annotated.pdf'),
                   ('DotPlot_annotation_markers', 'All_subgroups_annotated_marker_DotPlots.pdf')]:
    writer = PdfWriter()
    for group in groups:
        source = output / group / f'{stem}.pdf'
        reader = PdfReader(source)
        assert len(reader.pages) == 1, source
        writer.append(reader)
    destination = output / name
    temporary = destination.with_suffix('.pdf.partial')
    assert not temporary.exists(), temporary
    with temporary.open('wb') as fh:
        writer.write(fh)
    assert len(PdfReader(temporary).pages) == len(groups)
    temporary.replace(destination)
    print('PDF_COLLECTION_OK', name, len(groups), 'pages')
