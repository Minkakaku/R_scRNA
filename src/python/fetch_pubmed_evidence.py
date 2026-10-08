"""保存给定 PMID 的 NCBI 文献元数据及可获取的 PMC 正文，供注释追溯。"""
import argparse
import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--pmids', required=True)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=True)

def fetch(db, ids):
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urllib.parse.urlencode(
        {'db': db, 'id': ids, 'retmode': 'xml'})
    with urllib.request.urlopen(url, timeout=45) as response:
        return response.read(), url

xml, source = fetch('pubmed', args.pmids)
(args.out / 'pubmed_records.xml').write_bytes(xml)
records = []
for article in ET.fromstring(xml).findall('PubmedArticle'):
    ids = {node.get('IdType'): node.text for node in article.findall('./PubmedData/ArticleIdList/ArticleId')}
    record = dict(pmid=article.findtext('./MedlineCitation/PMID'),
                  title=''.join(article.find('.//ArticleTitle').itertext()),
                  journal=article.findtext('.//Journal/Title'),
                  year=article.findtext('.//JournalIssue/PubDate/Year'),
                  doi=ids.get('doi'), pmcid=ids.get('pmc'),
                  pubmed_url='https://pubmed.ncbi.nlm.nih.gov/' + ids['pubmed'] + '/',
                  abstract='\n'.join(''.join(node.itertext()) for node in article.findall('.//Abstract/AbstractText')),
                  retrieved_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), pubmed_api_url=source)
    if record['pmcid']:
        time.sleep(.4)
        try:
            full, full_url = fetch('pmc', record['pmcid'].replace('PMC', ''))
            root = ET.fromstring(full)
            if root.find('.//body') is None:
                raise ValueError('PMC response contains no full-text body')
            (args.out / (record['pmcid'] + '.xml')).write_bytes(full)
            paragraphs = [' '.join(''.join(node.itertext()).split()) for node in root.findall('.//body//p')]
            (args.out / (record['pmcid'] + '_paragraphs.txt')).write_text('\n\n'.join(paragraphs), encoding='utf-8')
            record['pmc_api_url'] = full_url
            record['full_text_saved'] = True
        except Exception as error:
            record['full_text_saved'] = False
            record['full_text_error'] = str(error)
    records.append(record)
    print(record['pmid'], record['title'], 'PMC_saved=' + str(record.get('full_text_saved')))
(args.out / 'references.json').write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding='utf-8')
assert {r['pmid'] for r in records} == set(args.pmids.split(','))
print('PUBMED_REFERENCE_RECORDS_OK', len(records))
