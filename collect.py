"""Public metadata ingestion. Standard library only; no patient-level data downloaded."""
import argparse
import concurrent.futures
import datetime as dt
import hashlib
import html
import json
import os
from pathlib import Path
import re
import threading
import time
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parent
TAXONOMY = json.loads((ROOT / 'taxonomy.json').read_text())
CURATION = json.loads((ROOT / 'curation.json').read_text()) if (ROOT / 'curation.json').exists() else {}
CACHE = ROOT / 'cache_cdx'
NOW = dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')
EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
lock = threading.Lock()
last_ncbi = 0.0
SPATIAL = r'\b(spatial (?:transcriptom|RNA|gene expression)|Visium|Xenium|CosMx|MERFISH|seqFISH|Slide.seq|Stereo.seq|GeoMx|digital spatial profiling|in situ sequencing)'
SCRNA = r'\b(scRNA[ -]?seq|single[ -]cell RNA[ -]?(?:seq|sequencing)|single[ -]cell sequencing)'
SNRNA = r'\b(snRNA[ -]?seq|single[ -](?:nucle\w+) (?:RNA|transcriptom))'
PLATFORMS = {'Xenium': ('Xenium', 'Imaging-based'), 'CosMx': ('CosMx', 'Imaging-based'), 'MERFISH': ('MERFISH', 'Imaging-based'), 'seqFISH': ('seqFISH', 'Imaging-based'), 'Visium HD': ('Visium HD', 'Sequencing-based'), 'Visium': ('Visium', 'Sequencing-based'), 'Slide-seq': ('Slide.seq', 'Sequencing-based'), 'Stereo-seq': ('Stereo.seq', 'Sequencing-based'), 'GeoMx': ('GeoMx|digital spatial profiling', 'Sequencing-based'), '10x Chromium': ('Chromium|10x (?:3|5)[\x27’]|10x Genomics single.cell', 'Sequencing-based'), 'Smart-seq': ('Smart.seq', 'Sequencing-based')}

def get(url, ttl=24*3600):
    global last_ncbi
    CACHE.mkdir(exist_ok=True)
    path = CACHE / (hashlib.sha256(url.encode()).hexdigest() + '.json')
    if path.exists() and time.time() - path.stat().st_mtime < ttl:
        return json.loads(path.read_text())
    for attempt in range(4):
        try:
            if 'ncbi.nlm.nih.gov' in url:
                with lock:
                    time.sleep(max(0, .38 - (time.monotonic() - last_ncbi)))
                    last_ncbi = time.monotonic()
            req = urllib.request.Request(url, headers={'User-Agent': 'HumanSpatialAtlas/1.0 (public research metadata; github.com/cemalyildiz/spatial_sc_atlas_cdx)'})
            with urllib.request.urlopen(req, timeout=45) as response:
                raw = response.read().decode('utf-8')
            value = json.loads(raw)
            if isinstance(value, dict) and ('error' in value or 'ERROR' in value):
                raise ValueError(str(value)[:180])
            path.write_text(json.dumps(value))
            return value
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)

def api(base, **params):
    return get(base + '?' + urllib.parse.urlencode(params))

def matches(text, term):
    return bool(re.search(r'(?<!\w)(?<!non-)(?<!non )' + re.escape(term) + r'(?!\w)', text, re.I))

def clean(text):
    return re.sub(r'\s+', ' ', html.unescape(re.sub('<[^>]+>', ' ', text or ''))).strip()

def classify(text):
    diseases, groups, tissues, subtypes = [], [], [], []
    for row in TAXONOMY:
        if any(matches(text, term) for term in row['terms']):
            diseases.append(row['disease']); groups.append(row['group']); tissues.append(row['tissue'])
            subtypes.extend(name for name, terms in row['subtypes'].items() if any(matches(text, t) for t in terms))
    if not diseases:
        if re.search(r'\b(cancer|carcinoma|tumou?r|malignan|neoplas)', text, re.I):
            diseases, groups = ['Other cancer'], ['Cancer']
        elif re.search(r'\b(disease|disorder|fibrosis|infection|injury|syndrome|inflamm)', text, re.I):
            diseases, groups = ['Other disease'], ['Other disease']
        else:
            return None
    mods = []
    if re.search(SPATIAL, text, re.I): mods.append('Spatial')
    if re.search(SCRNA, text, re.I): mods.append('scRNA-seq')
    if 'scRNA-seq' not in mods and (re.search(r'single[ -]cell and spatial (?:transcriptom|RNA)', text, re.I) or ('Spatial' not in mods and re.search(r'single[ -]cell transcriptom',text,re.I))): mods.append('scRNA-seq')
    if re.search(SNRNA, text, re.I): mods.append('snRNA-seq')
    platforms, approaches = [], []
    for name, (pattern, approach) in PLATFORMS.items():
        if re.search(pattern, text, re.I): platforms.append(name); approaches.append(approach)
    if 'Visium HD' in platforms and 'Visium' in platforms: platforms.remove('Visium')
    # Cell sequencing does not determine the spatial assay's approach.
    if 'Spatial' not in mods and any(m in mods for m in ['scRNA-seq', 'snRNA-seq']):
        approaches.append('Sequencing-based')
    return dict(diseases=sorted(set(diseases)), groups=sorted(set(groups)), tissues=sorted(set(tissues)) or ['Not specified'], subtypes=sorted(set(subtypes)) or ['Not specified'], modalities=mods, platforms=platforms or ['Not specified'], approaches=sorted(set(approaches)) or ['Not specified'])

def record(identifier, title, summary, url, source, released, **extra):
    if re.search(r'\[(?:bulk|ATAC|ChIP|WGS|WES)[^\]]*\]|\((?:bulk RNA|ATAC|ChIP)[^)]*\)',title,re.I): return None
    info = classify(title + ' ' + summary)
    if info and identifier in CURATION:
        info.update({k:v for k,v in CURATION[identifier].items() if k in info})
    if not info or not info['modalities']: return None
    both = 'Spatial' in info['modalities'] and len(info['modalities']) > 1
    text = title + '. ' + summary
    evidence = [s.strip()[:420] for s in re.split(r'(?<=[.!?])\s+', text) if re.search(SPATIAL+'|'+SCRNA+'|'+SNRNA, s, re.I)][:3]
    rec = dict(id=identifier, title=clean(title), summary=clean(summary)[:1300], sources=[source], url=url, released=(released or '')[:10].replace('/', '-'), updated='', checked=NOW, organism='Homo sapiens', pairing='Both reported — unverified' if both else 'Not applicable', validation='Metadata-derived', evidence=evidence, sample_count=None, cell_count=None, access='See source', availability=[], publications=[], related=[], source_records=[{'id':identifier,'source':source,'url':url}], **info)
    rec.update(extra)
    return rec

def geo_search(term, limit):
    result = api(EUTILS+'esearch.fcgi', db='gds', term=f'({term}) AND "Homo sapiens"[Organism] AND gse[Entry Type]', retmode='json', retmax=limit)['esearchresult']
    return result['idlist'], int(result['count'])

def geo_fetch(ids):
    records = []
    ids = list(dict.fromkeys(ids))
    for offset in range(0, len(ids), 60):
        raw = api(EUTILS+'esummary.fcgi', db='gds', id=','.join(ids[offset:offset+60]), retmode='json')['result']
        for uid in raw.get('uids', []):
            row = raw[uid]
            # Reject mixed-species series instead of silently labeling the entire study human.
            if row.get('taxon', '').strip() != 'Homo sapiens' or row.get('entrytype') != 'GSE': continue
            acc = row['accession']
            rec = record(acc, row['title'], row.get('summary',''), 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+acc, 'GEO', row.get('pdat',''), sample_count=int(row.get('n_samples') or len(row.get('samples',[]))) or None)
            if not rec: continue
            rec['publications'] = [{'pmid':str(p), 'url':'https://pubmed.ncbi.nlm.nih.gov/'+str(p)+'/'} for p in row.get('pubmedids',[])]
            if row.get('suppfile'): rec['availability'].append('Supplementary files')
            if any(x in row.get('suppfile','').upper() for x in ['H5','MTX','RDS','TXT','CSV']): rec['availability'].append('Processed files listed')
            rec['file_types'] = row.get('suppfile','')
            rec['sample_examples'] = [x.get('title','') for x in row.get('samples',[])[:5]]
            records.append(rec)
        print('GEO summaries:', min(offset+60,len(ids)), '/', len(ids), flush=True)
    return records

def collect_geo(limit):
    queries = []
    spatial = '"spatial transcriptomics" OR "spatial transcriptome" OR "spatial gene expression" OR Xenium OR Visium OR MERFISH OR CosMx OR "Slide-seq" OR "Stereo-seq" OR GeoMx'
    single = '"single cell RNA" OR "single-cell RNA" OR "single nucleus RNA" OR "single-nucleus RNA" OR "scRNA-seq" OR "snRNA-seq"'
    ids = ['200176078','200199102','200284230','200308146']
    for label, term, cap in [('Spatial methods',spatial,limit), ('Recent single-cell/nucleus',f'({single}) AND ("2024/01/01"[PDAT] : "3000"[PDAT])',limit)]:
        found,total = geo_search(term,cap); ids += found
        queries.append(dict(label=label,query=term,total=total,retrieved=len(found),limit=cap))
    # Independent subtype searches catch studies whose title omits the parent disease.
    for row in TAXONOMY:
        terms = list(dict.fromkeys(row['terms'] + [t for values in row['subtypes'].values() for t in values]))
        disease = ' OR '.join('"'+t+'"' for t in terms)
        term = f'({disease}) AND ({spatial} OR {single})'
        found,total = geo_search(term,35); ids += found
        queries.append(dict(label=row['disease']+' and subtypes',query=term,total=total,retrieved=len(found),limit=35))
    return geo_fetch(ids), dict(queries=queries, scanned=len(set(ids)), note='Spatial discovery plus recent single-cell studies and targeted disease/subtype searches. Query limits are shown; this is not an exhaustive census. Mixed-species series excluded.')

def walk(value):
    if isinstance(value, dict):
        yield value
        for child in value.values(): yield from walk(child)
    elif isinstance(value, list):
        for child in value: yield from walk(child)

def collect_arrayexpress(limit):
    queries, hits = [], {}
    for term in ['(spatial OR Visium OR Xenium OR CosMx OR MERFISH) AND "Homo sapiens"', '("single cell" OR "single-cell" OR "single nucleus") AND "Homo sapiens"']:
        scanned=0; total=0
        for page in range(1, max(2,limit//100+1)):
            response=api('https://www.ebi.ac.uk/biostudies/api/v1/ArrayExpress/search',query=term,pageSize=100,page=page,sortBy='release_date',sortOrder='descending')
            total=response['totalHits']
            for row in response.get('hits',[]): hits[row['accession']]=row
            scanned+=len(response.get('hits',[]))
            if scanned>=min(limit,total): break
        queries.append(dict(query=term,total=total,retrieved=scanned,limit=limit))
    records=[]; errors=[]
    def detail(hit):
        acc=hit['accession']
        if acc.startswith('E-GEOD-'): return None # canonical GEO record covers mirrored studies
        raw=get('https://www.ebi.ac.uk/biostudies/api/v1/studies/'+acc)
        attrs=list(walk(raw))
        organisms={a.get('value') for a in attrs if a.get('name','').lower()=='organism'}
        if organisms != {'Homo sapiens'}: return None
        section=raw.get('section',{})
        summary=' '.join(a.get('value','') for a in section.get('attributes',[]) if a.get('name') in ['Description','Study type'])
        rec=record(acc,hit['title'],summary,'https://www.ebi.ac.uk/biostudies/arrayexpress/studies/'+acc,'ArrayExpress',hit.get('release_date'))
        if rec:
            if hit.get('files'): rec['availability']=['Supplementary files']
            pmids=list(dict.fromkeys(a.get('value') for a in attrs if a.get('name','').lower() in ['pubmed id','pubmed'] and str(a.get('value','')).isdigit()))
            rec['publications']=[{'pmid':p,'url':'https://pubmed.ncbi.nlm.nih.gov/'+p+'/'} for p in pmids]
        return rec
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        futures={pool.submit(detail,h):h['accession'] for h in hits.values()}
        for future in concurrent.futures.as_completed(futures):
            try:
                item=future.result()
                if item: records.append(item)
            except Exception as e: errors.append(futures[future]+': '+str(e)[:120])
    return records, dict(queries=queries,scanned=len(hits),errors=errors,note='ArrayExpress collection only; GEO mirrors excluded. Organism verified from study metadata.')

def collect_cellxgene():
    collections=get('https://api.cellxgene.cziscience.com/curation/v1/collections')
    records=[]
    for col in collections:
        datasets=[d for d in col.get('datasets',[]) if {o['ontology_term_id'] for o in d.get('organism',[])}=={'NCBITaxon:9606'}]
        def dataset_modes(d):
            labels=[x['label'] for x in d.get('assay',[])]
            spatial=any(re.search(SPATIAL,a,re.I) for a in labels)
            rna=any(not re.search('ATAC|ChIP|protein|methyl',a,re.I) and re.search(r"10x|RNA|Smart.seq|Drop.seq|Seq.Well|sci.RNA|CEL.seq",a,re.I) and not re.search(SPATIAL,a,re.I) for a in labels)
            modes=['Spatial'] if spatial else []
            if rna and 'cell' in d.get('suspension_type',[]):modes.append('scRNA-seq')
            if rna and 'nucleus' in d.get('suspension_type',[]):modes.append('snRNA-seq')
            return modes
        diseased=[d for d in datasets if dataset_modes(d) and any(x.get('label') not in ['normal','unknown','na'] for x in d.get('disease',[]))]
        if not diseased: continue
        conditions=sorted({x['label'] for d in diseased for x in d.get('disease',[]) if x['label'] not in ['normal','unknown','na']})
        assays=sorted({x['label'] for d in diseased for x in d.get('assay',[])})
        susp={s for d in diseased for s in d.get('suspension_type',[])}
        mods=sorted({m for d in diseased for m in dataset_modes(d)})
        if not mods: continue
        title=col.get('name', 'CELLxGENE collection')
        summary=clean(col.get('description',''))
        info=classify(' '.join(conditions))
        if not info: info=classify(title+' '+summary)
        if not info: continue
        rec=record('CXG:'+col['collection_id'],title,summary+' '+ ' '.join(conditions)+' single-cell RNA sequencing',col['collection_url'],'CELLxGENE',col.get('published_at') or col.get('created_at'),modalities=mods,validation='Repository assay metadata',updated=(col.get('revised_at') or '')[:10])
        if not rec: continue
        for k in ['diseases','groups','subtypes']: rec[k]=info[k]
        rec['summary']=summary[:1300]
        rec['tissues']=sorted({x['label'] for d in diseased for x in d.get('tissue',[])})
        rec['platforms']=assays
        rec['approaches']=sorted({approach for name,(pattern,approach) in PLATFORMS.items() if any(re.search(pattern,a,re.I) for a in assays)}) or (['Sequencing-based'] if 'Spatial' not in mods else ['Not specified'])
        rec['availability']=['Processed matrix','Cell annotations']
        rec['access']='Open processed data'
        rec['pairing']='Same collection' if 'Spatial' in mods and len(mods)>1 else 'Not applicable'
        rec['evidence']=['Repository assays: '+', '.join(assays), 'Human dataset suspension types: '+', '.join(sorted(susp)), 'Repository disease annotations: '+', '.join(conditions)]
        rec['dataset_count']=len(diseased)
        rec['source_records']=[{'id':d['dataset_id'],'source':'CELLxGENE','url':col['collection_url']} for d in diseased]
        doi=col.get('doi')
        if doi: rec['publications']=[{'doi':doi,'url':'https://doi.org/'+doi.removeprefix('https://doi.org/')}]
        rec['links']=[{'label':l.get('link_name') or l.get('link_type','Source'),'url':l['link_url']} for l in col.get('links',[]) if l.get('link_url','').startswith('https://')]
        records.append(rec)
    return records,dict(scanned=len(collections),note='Human disease datasets only. Healthy control datasets may accompany the original collection. Collection membership does not establish donor matching.')

def collect_publications(limit):
    query='(TITLE_ABS:"spatial transcriptomics" OR TITLE_ABS:"single-cell RNA sequencing" OR TITLE_ABS:Xenium OR TITLE_ABS:CosMx) AND (TITLE_ABS:human OR TITLE_ABS:patients) AND FIRST_PDATE:[2024-01-01 TO 3000-12-31] sort_date:y'
    result=api('https://www.ebi.ac.uk/europepmc/webservices/rest/search',query=query,format='json',pageSize=min(limit,1000),resultType='core')
    pubs=[]; geo=set()
    for p in result.get('resultList',{}).get('result',[]):
        abstract=clean(p.get('abstractText','')); text=p.get('title','')+' '+abstract
        if not classify(text): continue
        accs=sorted(set(re.findall(r'\bGSE\d+\b',text)))
        geo.update(accs)
        pubs.append(dict(pmid=p.get('pmid'),doi=p.get('doi'),title=p.get('title'),date=p.get('firstPublicationDate'),url='https://europepmc.org/article/'+p.get('source','MED')+'/'+p['id'],accessions=accs))
    return pubs,geo,dict(scanned=len(result.get('resultList',{}).get('result',[])),total=result.get('hitCount'),query=query,note='Publication discovery is separate from dataset availability. Papers without verified repository records are not counted as datasets.')

def link_publications(records,pubs):
    by_pmid={p['pmid']:p for p in pubs if p.get('pmid')}
    by_doi={p['doi'].lower():p for p in pubs if p.get('doi')}
    groups={}
    for rec in records:
        for pub in rec['publications']:
            extra=by_pmid.get(pub.get('pmid')) or by_doi.get(pub.get('doi','').lower())
            if extra: pub.update(extra)
            key=('pmid:'+pub['pmid']) if pub.get('pmid') else ('doi:'+pub.get('doi','').lower())
            if key not in ['doi:','pmid:']: groups.setdefault(key,[]).append(rec)
    for siblings in groups.values():
        for rec in siblings:
            rec['related']=list(dict.fromkeys(rec['related']+[s['id'] for s in siblings if s['id']!=rec['id']]))
    # Only explicit accession links are merged, never a shared title or PMID alone.
    geo_map={r['id']:r for r in records if r['id'].startswith('GSE')}
    for rec in records:
        for link in rec.get('links',[]):
            for acc in re.findall(r'\bGSE\d+\b',link['url']):
                if acc in geo_map:
                    target=geo_map[acc]
                    rec['related']=list(dict.fromkeys(rec['related']+[acc]))
                    target['related']=list(dict.fromkeys(target['related']+[rec['id']]))
    return records

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--limit',type=int,default=1500); parser.add_argument('--ae-limit',type=int,default=200)
    args=parser.parse_args()
    old=json.loads((ROOT/'catalog.json').read_text()) if (ROOT/'catalog.json').exists() else {'records':[], 'sources':{}}
    records=[]; statuses={}; pubs=old.get('publications',[])
    tasks={'GEO':lambda:collect_geo(args.limit),'ArrayExpress':lambda:collect_arrayexpress(args.ae_limit),'CELLxGENE':collect_cellxgene}
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        futures={pool.submit(fn):name for name,fn in tasks.items()}
        for future in concurrent.futures.as_completed(futures):
            name=futures[future]
            try:
                rows,coverage=future.result()
                if not rows: raise ValueError('Empty source result; retaining previous records')
                prior={r['id']:r for r in old['records'] if name in r['sources'] and (name=='CELLxGENE' or record(r['id'],r['title'],r['summary'],r['url'],name,r['released']))}
                prior.update({r['id']:r for r in rows})
                records.extend(prior.values())
                statuses[name]=dict(status='partial' if coverage.get('errors') else 'ok',last_success=NOW,count=len(prior),coverage=coverage)
                print(name, len(rows), 'records',flush=True)
            except Exception as e:
                rows=[r for r in old['records'] if name in r['sources']]; records.extend(rows)
                statuses[name]=dict(status='failed',last_success=old.get('sources',{}).get(name,{}).get('last_success'),count=len(rows),error=str(e)[:240])
                print(name,'FAILED',str(e),flush=True)
    try:
        pubs,discovered,coverage=collect_publications(300)
        known={r['id'] for r in records}
        if discovered-known: records.extend(geo_fetch(['200'+acc[3:] for acc in discovered-known]))
        statuses['PubMed / Europe PMC']=dict(status='ok',last_success=NOW,count=len(pubs),coverage=coverage)
    except Exception as e:
        statuses['PubMed / Europe PMC']=dict(status='failed',count=len(pubs),last_success=old.get('sources',{}).get('PubMed / Europe PMC',{}).get('last_success'),error=str(e)[:240])
    records=list({r['id']:r for r in records}.values())
    records=link_publications(records,pubs)
    curated=CURATION
    for rec in records:
        if rec['id'] in curated: rec.update(curated[rec['id']])
    for name, status in statuses.items():
        if name != 'PubMed / Europe PMC': status['count']=sum(name in r['sources'] for r in records)
    records.sort(key=lambda r:r['released'],reverse=True)
    payload=dict(schema_version=1,generated_at=NOW,sources=statuses,records=records,publications=pubs,taxonomy=TAXONOMY)
    if not records: raise RuntimeError('No records; existing catalogue untouched')
    temp=ROOT/'catalog.tmp'; temp.write_text(json.dumps(payload,ensure_ascii=False,separators=(',',':'))); temp.replace(ROOT/'catalog.json')
    print('Saved',len(records),'records',flush=True)

if __name__=='__main__': main()
