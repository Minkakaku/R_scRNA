"""核验评分和组成交付文件，保护已验证的原结果，并渲染合辑用于目视检查。"""
import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

from PIL import Image, ImageDraw
from pypdf import PdfReader


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest().lower()


def rows(path, delimiter='\t'):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def main(project, config_path):
    project = project.resolve()
    config = json.loads(config_path.read_text(encoding='utf-8'))
    scoring = config['neutrophil_scoring_request']
    composition = config['composition_request']
    assert scoring['status'] == composition['status'] == 'completed'
    assert scoring['run_scoring'] and scoring['confirmed_on']
    assert sha(project / scoring['input']) == scoring['input_sha256']
    assert sha(project / scoring['mapping']) == scoring['mapping_sha256']
    assert sha(project / scoring['user_reference']) == scoring['reference_sha256']
    record = rows(project / scoring['evidence_dir'] / 'object_validation.tsv')[0]
    assert sha(project / scoring['output']) == record['output_sha256']
    numeric = json.loads((project / scoring['evidence_dir'] / 'independent_validation.json').read_text())
    counts = json.loads((project / composition['evidence_dir'] / 'independent_validation.json').read_text())
    roe = json.loads((project / composition['evidence_dir'] / 'roe_dotplot_validation.json').read_text())
    assert numeric['maximum_score_error'] < 1e-10 and counts['all_counts_fractions_ROE_verified']
    assert roe['zero_ROE_has_no_point'] and roe['ROE_values_unchanged']
    assert roe['max_formula_error'] < 1e-12
    archived = protected = 0
    for name in ('73_scoring_archive.csv', '77_composition_archive.csv', '80_roe_heatmap_archive.csv', '82_roe_dotplot_archive.csv', '85_validation_code_archive.csv'):
        for row in rows(project / 'logs' / name, ','):
            assert sha(row['archive']) == row['sha256'].lower(), row['archive']
            archived += 1
    for name in ('73_scoring_protected_results.csv', '80_roe_protected_results.csv'):
        for row in rows(project / 'logs' / name, ','):
            assert sha(row['path']) == row['sha256'].lower(), row['path']
            protected += 1
    results = project / 'results'
    pdfs = sorted(results.rglob('*.pdf'))
    pngs = sorted(results.rglob('*.png'))
    assert len(list(results.rglob('README*'))) == 1
    for path in pdfs:
        check = subprocess.run(['pdfimages', '-list', str(path)], capture_output=True, text=True, check=True)
        assert not re.search(r'^\s*\d+\s+\d+\s+(?:image|mask|smask)\b', check.stdout, flags=re.M), path
        if len(PdfReader(path).pages) == 1:
            assert path.with_suffix('.png').exists(), path
    for path in pngs:
        with Image.open(path) as im:
            assert all(abs(v - 600) < 1 for v in im.info.get('dpi', (0, 0))), path
        assert path.with_suffix('.pdf').exists()
    readme = (results / 'README.md').read_text(encoding='utf-8')
    for target in re.findall(r'\]\(([^)]+)\)', readme):
        if not re.match(r'^[a-z]+://', target) and not target.startswith('#'):
            assert (results / target.split('#')[0]).exists(), target
    scopes = ['Main_lineages'] + composition['scopes']
    new_pdfs = []
    for scope in scopes:
        folder = project / (composition['main_results'] if scope == 'Main_lineages' else f"{composition['subgroup_results']}/{scope}")
        assert (folder / 'Composition_ROE.tsv').exists()
        for kind in ('sample', 'group'):
            for plot in ('Counts', 'Proportions', 'ROE'):
                path = folder / f'{plot}_by_{kind}.pdf'
                assert path.exists() and path.with_suffix('.png').exists()
                new_pdfs.append(path)
    for name in ('Heatmap_figure_scores', 'Scores_figure_by_sample', 'DotPlot_figure_score_genes', 'UMAP_figure_scores'):
        new_pdfs.append(project / scoring['results_dir'] / f'{name}.pdf')
    for path in new_pdfs:
        text = '\n'.join(page.extract_text() for page in PdfReader(path).pages)
        assert not re.search(r'\bC\d+\b', text), path
    book = project / composition['main_results'] / 'All_composition_ROE.pdf'
    assert len(PdfReader(book).pages) == 2 * len(scopes)
    qa = project / 'logs/85_figure_qa'
    qa.mkdir(exist_ok=True)
    subprocess.run(['pdftoppm', '-r', '85', '-png', str(book), str(qa / 'composition')], check=True, capture_output=True)
    pages = sorted(qa.glob('composition-*.png'))
    assert len(pages) == 2 * len(scopes)
    for start in range(0, len(pages), 8):
        sheet = Image.new('RGB', (1600, 1200), 'white')
        draw = ImageDraw.Draw(sheet)
        for j, path in enumerate(pages[start:start + 8]):
            with Image.open(path) as im:
                im.thumbnail((790, 270))
                x, y = (j % 2) * 800, (j // 2) * 300
                sheet.paste(im, (x, y + 20))
                draw.text((x + 5, y + 3), f'Page {start + j + 1}', fill='black')
        sheet.save(qa / f'contact_{start // 8 + 1}.png')
    output = {
        'score_cells': numeric['cells'], 'score_clusters': numeric['clusters'],
        'composition_cells': counts['cells'], 'composition_plot_pairs': 42,
        'score_plot_pairs': 4, 'ROE_dotplots': roe['figures'], 'overview_pages': len(pages),
        'all_formal_PDFs': len(pdfs), 'all_formal_PNGs': len(pngs),
        'all_PDFs_vector': True, 'all_PNGs_600dpi': True,
        'archived_files_verified': archived, 'protected_checks_verified': protected,
        'input_scored_QS_and_annotation_hashes_verified': True,
        'README_links_and_new_figure_labels_verified': True,
        'new_figure_hashes': {str(p.relative_to(project)): sha(p) for p in new_pdfs + [book]},
    }
    (project / 'logs/85_delivery_validation.json').write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in output.items() if k != 'new_figure_hashes'}, ensure_ascii=False, indent=2))
    print('SCORES_COMPOSITION_DELIVERY_VALIDATION_OK')


if __name__ == '__main__':
    main(Path(sys.argv[1]), Path(sys.argv[2]))
