"""Validate the published catalogue contract without requiring network access."""
import datetime as dt
import json
from pathlib import Path
from urllib.parse import urlparse

def validate():
    data=json.loads((Path(__file__).parent/'catalog.json').read_text())
    rows=data['records']
    assert len(rows)>0, 'No dataset records'
    assert len({r['id'] for r in rows})==len(rows), 'Duplicate accessions'
    for r in rows:
        assert r['organism']=='Homo sapiens', r['id']
        assert r['title'] and r['sources'] and r['modalities'], r['id']
        assert set(r['modalities']) <= {'Spatial','scRNA-seq','snRNA-seq'}, r['id']
        assert urlparse(r['url']).scheme=='https', r['id']
        for field in ['groups','diseases','subtypes','tissues','approaches','platforms','evidence','publications','related','availability']:
            assert isinstance(r[field],list),(r['id'],field)
        if r['released']:
            assert dt.date.fromisoformat(r['released']) <= dt.datetime.now(dt.timezone.utc).date(), r['id']
        if r['pairing'] in ['Matched donors / samples','Same study — verified']:
            assert r['validation']=='Source-reviewed pairing' and r.get('evidence_url') and r['evidence'], r['id']
        if r['pairing']!='Not applicable':
            assert 'Spatial' in r['modalities'] and len(r['modalities'])>1,r['id']
    for source in ['GEO','ArrayExpress','CELLxGENE','PubMed / Europe PMC']:
        assert source in data['sources'],source
    print(f'Validated {len(rows):,} human records, unique accessions, dates, links and pairing evidence.')

if __name__=='__main__': validate()
