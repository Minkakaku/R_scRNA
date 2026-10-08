"""按既有哈希清单恢复正式结果；先归档被替换文件，保留完整历史。"""
import csv
import hashlib
import shutil
import sys
from pathlib import Path

project = Path(sys.argv[1]).resolve()
manifest = project / sys.argv[2]
archive = project / sys.argv[3]
report = project / sys.argv[4]

def sha256(path):
    with path.open('rb') as fh:
        return hashlib.file_digest(fh, 'sha256').hexdigest().upper()

with manifest.open(encoding='utf-8-sig', newline='') as fh:
    rows = list(csv.DictReader(fh))
assert rows and not archive.exists(), '需要非空清单和新的归档目录'
for row in rows:
    target = (project / 'results' / row['path']).resolve()
    assert target.is_relative_to(project / 'results') and target.is_file()
    assert sha256(Path(row['archive'])) == row['sha256'].upper()

# 先完整归档本轮将替换的结果及说明，再恢复选定版本。
paths = [Path('results') / r['path'] for r in rows] + [Path(p) for p in sys.argv[5:]]
saved = []
for relative in paths:
    source = (project / relative).resolve()
    destination = archive / relative
    assert source.is_relative_to(project) and source.is_file()
    destination.parent.mkdir(parents=True, exist_ok=True)
    digest = sha256(source)
    shutil.copy2(source, destination)
    assert sha256(destination) == digest
    saved.append({'path': str(relative), 'sha256': digest, 'archive': str(destination)})
with report.open('w', encoding='utf-8', newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=['path', 'sha256', 'archive'])
    writer.writeheader()
    writer.writerows(saved)
for row in rows:
    target = project / 'results' / row['path']
    shutil.copy2(row['archive'], target)
    assert sha256(target) == row['sha256'].upper()
print(f'RESULTS_RESTORED_OK restored={len(rows)} archived={len(saved)}')
