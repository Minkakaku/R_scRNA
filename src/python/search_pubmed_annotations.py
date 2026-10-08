"""按项目检索式保存PubMed检索结果和候选论文摘要，保留零命中与失败记录。"""
import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--queries', required=True, type=Path)
parser.add_argument('--out', required=True, type=Path)
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=False)


def get(endpoint, params):
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/' + endpoint + '?' + urllib.parse.urlencode(params)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                data = response.read()
            time.sleep(.4)
            return data, url
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))


queries = list(csv.DictReader(args.queries.open(encoding='utf-8-sig'), delimiter='\t'))
searches, pmids = [], set()
for row in queries:
    raw, url = get('esearch.fcgi', dict(db='pubmed', term=row['query'], retmode='json', retmax=6, sort='relevance'))
    data = json.loads(raw)['esearchresult']
    searches.append(dict(row, source_url=url, retrieved_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), result=data))
    pmids.update(data['idlist'])
    (args.out / 'searches.json').write_text(json.dumps(searches, ensure_ascii=False, indent=2), encoding='utf-8')
    print(row['topic'], 'hits=' + data['count'], ','.join(data['idlist']), flush=True)
records = []
ids = sorted(pmids)
for start in range(0, len(ids), 40):
    raw, url = get('efetch.fcgi', dict(db='pubmed', id=','.join(ids[start:start+40]), retmode='xml'))
    (args.out / f'pubmed_candidates_{start}.xml').write_bytes(raw)
    for article in ET.fromstring(raw).findall('PubmedArticle'):
        identifiers = {n.get('IdType'): n.text for n in article.findall('./PubmedData/ArticleIdList/ArticleId')}
        pmid = article.findtext('./MedlineCitation/PMID')
        records.append(dict(pmid=pmid, title=''.join(article.find('.//ArticleTitle').itertext()),
            abstract='\n'.join(''.join(n.itertext()) for n in article.findall('.//Abstract/AbstractText')),
            journal=article.findtext('.//Journal/Title'), year=article.findtext('.//JournalIssue/PubDate/Year'),
            publication_types=[n.text for n in article.findall('.//PublicationTypeList/PublicationType')],
            doi=identifiers.get('doi'), pmcid=identifiers.get('pmc'),
            pubmed_url=f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/', pubmed_api_url=url))
(args.out / 'references.json').write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding='utf-8')
assert set(ids) == {r['pmid'] for r in records}
print('PUBMED_SEARCH_COMPLETE', len(queries), 'queries', len(records), 'records', flush=True)
