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
    from assay_scope import scoped
    from collect import CURATION
    for r in rows:
        expected=scoped(r,CURATION.get(r['id'],{}).get('accession_assay_review'))
        for k,v in expected.items():assert r.get(k)==v,(r['id'],k,'unsupported accession label')
        assert r['organism']=='Homo sapiens', r['id']
        assert r['title'] and r['sources'] and r['modalities'], r['id']
        assert set(r['modalities']) <= {'Spatial','scRNA-seq','snRNA-seq'}, r['id']
        assert urlparse(r['url']).scheme=='https', r['id']
        for field in ['groups','diseases','subtypes','tissues','approaches','platforms','evidence','publications','related','availability']:
            assert isinstance(r[field],list),(r['id'],field)
        obs=r['sample_observations']
        assert set(s['id'] for s in obs)<=set(s['id'] for s in r.get('sample_metadata',[])),r['id']
        if 'CELLxGENE' not in r['sources']:
            for field,key in [('sampled_tissues','tissues'),('sampled_diagnoses','diagnoses'),('clinical_groups','clinical'),('sample_cell_types','cell_types')]:
                assert set(r[field])=={v for s in obs for v in s[key]},(r['id'],field)
        assert set(r['file_categories'])==set(r['file_inventory']),r['id']
        for files in r['file_inventory'].values():
            assert all(urlparse(f['url']).scheme=='https' for f in files),r['id']
        if r['released']:
            assert dt.date.fromisoformat(r['released']) <= dt.datetime.now(dt.timezone.utc).date(), r['id']
        if r['pairing'] in ['Matched donors / samples','Same study — verified']:
            assert r['validation'] in ['Source-reviewed pairing','Repository relationship evidence'] and r.get('evidence_url') and r['evidence'], r['id']
            if r['pairing']=='Matched donors / samples':assert r['validation']=='Source-reviewed pairing',r['id']
        if r['pairing']!='Not applicable':
            assert 'Spatial' in r['modalities'] and len(r['modalities'])>1,r['id']
    ids={r['id'] for r in rows}
    for edge in data.get('relationships',[]):
        assert len(edge['accessions'])==2 and all(a in ids for a in edge['accessions']),edge['id']
        if edge['status']=='Matched donors / samples':assert edge['reviewed'] and edge['matched_sample_count'],edge['id']
    for source in ['GEO','ArrayExpress','CELLxGENE','PubMed / Europe PMC']:
        assert source in data['sources'],source
    for publication in data.get('publications',[]):
        if publication.get('date'):
            assert dt.date.fromisoformat(publication['date']) <= dt.datetime.now(dt.timezone.utc).date(), publication.get('pmid')
    print(f'Validated {len(rows):,} human records, unique accessions, dates, links and pairing evidence.')

if __name__=='__main__': validate()
