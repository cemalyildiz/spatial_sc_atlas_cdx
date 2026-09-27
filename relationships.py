"""Evidence-preserving accession relationships; public metadata, never patient files."""
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations
import hashlib
import json
import re
import time
import urllib.request

METHOD = r'(?:scRNA[ -]?seq|snRNA[ -]?seq|single[ -]cell RNA[ -]?seq(?:uencing)?|single[ -]nucleus RNA[ -]?seq(?:uencing)?|spatial transcriptomics|spatial|(?<![A-Za-z])ST(?![A-Za-z])|Visium|Xenium|CosMx)'

def components(title, samples=()):
    """Deposited assay evidence from the component title, then sample labels, not background."""
    import collect as c
    def detect(text):
        mods=[]
        if re.search(c.SPATIAL, text, re.I) or re.search(r'\[ST\]|\(ST\)',text): mods.append('Spatial')
        if re.search(c.SCRNA,text,re.I): mods.append('scRNA-seq')
        if re.search(c.SNRNA,text,re.I): mods.append('snRNA-seq')
        return mods
    brackets=re.findall(r'\[([^]]+)\]',title)
    explicit=detect(' '.join('['+s+']' for s in brackets))
    return explicit or detect(title) or sorted({m for s in samples for m in detect(s.get('title',''))})

def title_key(title):
    title=re.sub(r'\[[^]]*\]|\([^)]*\)', ' ',title)
    title=re.sub(METHOD,' ',title,flags=re.I)
    return re.sub(r'\W+',' ',title).strip().lower()

def sample_key(title):
    s=re.sub(METHOD,' ',title,flags=re.I)
    s=re.sub(r'[^A-Za-z0-9]+','',s).upper()
    # Generic labels and pure numbers are not donor identifiers.
    return s if re.search('[A-Z]',s) and re.search('[0-9]',s) and len(s)>=3 else None

def add_geo_metadata(rec,row):
    samples=[{'id':s.get('accession',''),'title':s.get('title','')} for s in row.get('samples',[])[:500]]
    rec['sample_metadata']=samples
    rec['sample_metadata_complete']=len(samples)==rec.get('sample_count')
    rec['accession_modalities']=components(rec['title'],samples)
    if rec['accession_modalities']:
        rec['modalities']=rec['accession_modalities'][:]
        rec['pairing']='Both reported — unverified' if 'Spatial' in rec['modalities'] and len(rec['modalities'])>1 else 'Not applicable'
    rec['bioprojects']=re.findall(r'PRJNA\d+',str(row.get('bioproject','')))


def series_detail(acc):
    import collect as c
    url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief'
    c.CACHE.mkdir(exist_ok=True)
    path=c.CACHE/('series-'+acc+'.json')
    if path.exists() and time.time()-path.stat().st_mtime<7*86400:return json.loads(path.read_text())
    for attempt in range(3):
        try:
            with c.lock:
                time.sleep(max(0,.38-(time.monotonic()-c.last_ncbi)));c.last_ncbi=time.monotonic()
            req=urllib.request.Request(url,headers={'User-Agent':'HumanSpatialAtlas/1.0 public metadata'})
            with urllib.request.urlopen(req,timeout=25) as response:text=response.read().decode('utf-8')
            if '^SERIES = '+acc not in text:raise ValueError('GEO series response unavailable')
            values=defaultdict(list)
            for line in text.splitlines():
                match=re.match(r'!Series_(overall_design|relation|pubmed_id|summary) = (.*)',line)
                if match:values[match[1]].append(match[2])
            # Persist only the permitted study fields, not submitter contact information.
            result={'overall_design':' '.join(values['overall_design']), 'repository_relations':values['relation'], 'design_checked':c.NOW}
            path.write_text(json.dumps(result));return result
        except Exception:
            if attempt==2:raise
            time.sleep(2**attempt)


def publication_keys(rec):
    keys=set()
    for p in rec.get('publications',[]):
        if p.get('pmid'):keys.add('PMID:'+str(p['pmid']))
        if p.get('doi'):keys.add('DOI:'+re.sub(r'^https?://(?:dx\.)?doi.org/','',p['doi']).lower())
    return keys


def resolve_publication_ids(records,stats):
    import collect as c
    pmids=sorted({str(p['pmid']) for r in records for p in r.get('publications',[]) if str(p.get('pmid','')).isdigit() and not p.get('doi')})
    doi_by_pmid={}
    for i in range(0,len(pmids),100):
        try:
            result=c.api(c.EUTILS+'esummary.fcgi',db='pubmed',id=','.join(pmids[i:i+100]),retmode='json')['result']
            for uid in result.get('uids',[]):
                for ident in result[uid].get('articleids',[]):
                    if ident.get('idtype')=='doi':doi_by_pmid[uid]=ident['value']
        except Exception as e:stats['errors'].append('PMID/DOI resolution: '+str(e)[:120])
    for r in records:
        for p in r.get('publications',[]):
            if not p.get('doi') and str(p.get('pmid')) in doi_by_pmid:p['doi']=doi_by_pmid[str(p['pmid'])]
    stats['publication_dois_resolved']=len(doi_by_pmid)


def pair_candidates(records):
    indexes=defaultdict(set);known={r['id'] for r in records};pairs=defaultdict(set)
    for r in records:
        for k in publication_keys(r):indexes[k].add(r['id'])
        for k in r.get('bioprojects',[]):indexes['BioProject:'+k].add(r['id'])
        for relation in r.get('repository_relations',[]):
            for k in re.findall(r'PRJNA\d+|SRP\d+|GSE\d+',relation):
                indexes['Repository:'+k].add(r['id'])
                if k in known and k!=r['id']:pairs[tuple(sorted([r['id'],k]))].add('Repository:'+k)
        for k in re.findall(r'\bGSE\d+\b',r.get('overall_design','')+' '+r.get('summary','')):
            if k in known and k!=r['id']:pairs[tuple(sorted([r['id'],k]))].add('Explicit accession mention')
        for link in r.get('links',[]):
            for k in re.findall(r'\bGSE\d+\b',link.get('url','')):
                if k in known and k!=r['id']:pairs[tuple(sorted([r['id'],k]))].add('Explicit source link')
    for key,ids in indexes.items():
        # Relations are pairwise; no transitive union that merges distinct studies.
        for a,b in combinations(sorted(ids),2):pairs[(a,b)].add(key)
    return pairs


def is_complementary(a,b):
    ma=set(a.get('accession_modalities',[]));mb=set(b.get('accession_modalities',[]))
    return ('Spatial' in ma and bool(mb&{'scRNA-seq','snRNA-seq'})) or ('Spatial' in mb and bool(ma&{'scRNA-seq','snRNA-seq'}))


def external_reference(a,b):
    if 'Spatial' not in a.get('accession_modalities',[]):return ''
    for sentence in re.split(r'(?<=[.!?])\s+',a.get('overall_design','')+' '+a.get('summary','')):
        if b['id'] in sentence and re.search(r'\breference\b',sentence,re.I) and re.search(r'\b(public(?:ly)?|published|external|downloaded|previous(?:ly)?)\b',sentence,re.I) and re.search(r'single.cell|scRNA|single.nucleus|snRNA',sentence,re.I):
            return sentence[:900]
    return ''


def build_relationships(records,curation):
    by_id={r['id']:r for r in records};edges=[]
    for r in records:
        r['relationship_ids']=[]
        r['relationship_levels']=[]
        if r.get('relationship_auto'):
            r.update(r.pop('relationship_base'))
            r.pop('relationship_auto',None)
    for (aid,bid),basis in sorted(pair_candidates(records).items()):
        a,b=by_id[aid],by_id[bid]
        curated=next((r for r in [a,b] if r['id'] in curation and (bid if r is a else aid) in curation[r['id']].get('related',[])),None)
        if not is_complementary(a,b) and not curated:continue
        same_title=title_key(a['title'])==title_key(b['title']) and len(title_key(a['title']))>35
        shared_publication=publication_keys(a)&publication_keys(b)
        same_design=bool(a.get('overall_design')) and a.get('overall_design')==b.get('overall_design') and len(a['overall_design'])>60
        shared_samples=sorted(set(filter(None,(sample_key(s['title']) for s in a.get('sample_metadata',[])))) & set(filter(None,(sample_key(s['title']) for s in b.get('sample_metadata',[])))))
        status='Related study — unverified';evidence=['Complementary accession-level assay labels; shared identifiers alone do not prove a common cohort.']
        ref=external_reference(a,b) or external_reference(b,a)
        if ref:
            status='External reference';evidence=[ref]
        elif same_title and same_design and shared_publication:
            status='Same study — verified';evidence=['The component records share a study title, an identical experimental design and a linked publication. Donor/sample matching has not been independently verified.',a['overall_design'][:1000]]
        if curated:
            status=curated['pairing'];evidence=[t.replace('This accession',curated['id']) for t in curated['evidence']]
            shared_samples=curated.get('paired_sample_ids',shared_samples)
        # Sample label overlap is displayed only as a candidate unless source reviewed.
        edge={'id':aid+'--'+bid,'accessions':[aid,bid],'status':status,'basis':sorted(basis),'evidence':evidence,'shared_sample_labels':shared_samples,'matched_sample_count':curated.get('paired_sample_count') if curated else None,'source_urls':[a['url'],b['url']], 'reviewed':bool(curated),'modalities':sorted(set(a.get('accession_modalities',[])+b.get('accession_modalities',[])))}
        edges.append(edge)
        for rec,other in [(a,b),(b,a)]:
            rec['relationship_ids'].append(edge['id'])
            rec['relationship_levels']=sorted(set(rec['relationship_levels']+[status]))
            rec['related']=sorted(set(rec.get('related',[])+[other['id']]))
            if status not in ['Same study — verified','External reference'] or rec['id'] in curation:continue
            # Leave accession_modalities intact; union is explicitly study-level.
            ranks={'Same study — verified':3,'External reference':2,'Both reported — unverified':1,'Not applicable':0}
            if ranks.get(rec['pairing'],4)>ranks[status]:continue
            if not rec.get('relationship_auto'):
                rec['relationship_base']={k:rec[k] for k in ['modalities','pairing','validation','evidence']}
                rec['relationship_auto']=True
            rec['modalities']=sorted(set(rec['modalities']+edge['modalities']))
            rec['pairing']=status;rec['validation']='Repository relationship evidence';rec['evidence']=evidence
            rec['evidence_url']=rec['url']
    return edges


def enrich(records,pubs,discover=True):
    import collect as c
    stats={'checked_at':c.NOW,'input_records':len(records),'errors':[],'discovery_limit':2500}
    known={r['id'] for r in records};found=set()
    if discover:
        for r in records:
            text=r.get('summary','')+' '+r.get('overall_design','')+' '+' '.join(l.get('url','') for l in r.get('links',[]))+' '+' '.join(r.get('repository_relations',[]))
            found.update('200'+acc[3:] for acc in re.findall(r'\bGSE\d+\b',text) if acc not in known)
    pmids=sorted({str(p['pmid']) for r in records for p in r.get('publications',[]) if str(p.get('pmid','')).isdigit()})
    if discover:
        for i in range(0,len(pmids),100):
            try:
                response=c.api(c.EUTILS+'elink.fcgi',dbfrom='pubmed',db='gds',id=','.join(pmids[i:i+100]),retmode='json')
                for group in response.get('linksets',[]):
                    for links in group.get('linksetdbs',[]):
                        found.update(uid for uid in links.get('links',[]) if uid.startswith('200') and 'GSE'+uid[3:] not in known)
            except Exception as e:stats['errors'].append('Publication discovery: '+str(e)[:140])
            print('Publication links:',min(i+100,len(pmids)),'/',len(pmids),flush=True)
    stats.update(publications_searched=len(pmids) if discover else 0,discovered_geo_accessions=len(found))
    # Refresh all indexed GEO summaries to recover assay/sample metadata consistently.
    ids=['200'+r['id'][3:] for r in records if re.fullmatch(r'GSE\d+',r['id'])]
    ids+=sorted(found)[:stats['discovery_limit']]
    refreshed={}
    for i in range(0,len(ids),60):
        try:
            for r in c.geo_fetch(ids[i:i+60]):refreshed[r['id']]=r
        except Exception as e:stats['errors'].append('GEO summary batch: '+str(e)[:140])
    merged={r['id']:r for r in records}
    # Retain existing design fields if fresh summaries do not carry them.
    for aid,r in refreshed.items():
        merged[aid]={**merged.get(aid,{}),**r}
        merged[aid].pop('relationship_auto',None)
        merged[aid].pop('relationship_base',None)
    records=list(merged.values());stats['new_records']=len(set(merged)-known)
    records=c.link_publications(records,pubs)
    resolve_publication_ids(records,stats)
    for r in records:
        if r['id'] in c.CURATION:r.update(c.CURATION[r['id']])
        if 'CELLxGENE' in r['sources']:r['accession_modalities']=r['modalities'][:]
        elif 'accession_modalities' not in r:r['accession_modalities']=components(r['title'],r.get('sample_metadata',[]))
    candidates=pair_candidates(records)
    # Read experimental design for linked accession candidates and records reporting both.
    target={aid for pair in candidates for aid in pair if aid.startswith('GSE')}
    target.update(r['id'] for r in records if r['id'].startswith('GSE') and 'Spatial' in r['modalities'] and len(r['modalities'])>1)
    by_id={r['id']:r for r in records};success=0
    pending=[]
    for aid in sorted(target):
        r=by_id[aid]
        checked=r.get('design_checked')
        if checked and r.get('overall_design') and (c.dt.datetime.now(c.dt.timezone.utc)-c.dt.datetime.fromisoformat(checked)).total_seconds()<7*86400:
            success+=1
        else:pending.append(aid)
    stats['designs_cached']=success
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(series_detail,aid):aid for aid in pending}
        for f in as_completed(futures):
            aid=futures[f]
            try:
                by_id[aid].update(f.result());success+=1
            except Exception as e:stats['errors'].append(aid+': '+str(e)[:120])
            if (success+len(stats['errors']))%50==0:print('Experimental designs:',success,'/',len(target),flush=True)
    for r in records:
        info=c.classify(r['title']+' '+r.get('summary','')+' '+r.get('overall_design',''))
        if info and r['id'] not in c.CURATION:
            for field in ['platforms','approaches']:
                values=sorted((set(r[field])|set(info[field]))-{'Not specified'})
                r[field]=values or ['Not specified']
    stats.update(designs_requested=len(target),designs_read=success,records_scanned=len(records))
    edges=build_relationships(records,c.CURATION)
    stats['relationship_counts']=dict(Counter(e['status'] for e in edges))
    stats['status']='partial' if stats['errors'] or len(found)>stats['discovery_limit'] else 'ok'
    stats['note']='Every indexed record was checked for direct publication, BioProject/SRA and accession links. Experimental designs were read for linked GEO candidates and records reporting both modalities. Shared sample labels are candidates, not proof of matched donors. Discovery is one hop, capped per run; unread or unavailable metadata remains unverified.'
    return records,edges,stats


def main():
    import collect as c
    path=c.ROOT/'catalog.json';data=json.loads(path.read_text())
    data['records'],data['relationships'],data['relationship_scan']=enrich(data['records'],data.get('publications',[]))
    data['records'].sort(key=lambda r:r['released'],reverse=True)
    for name,status in data['sources'].items():
        if name!='PubMed / Europe PMC':status['count']=sum(name in r['sources'] for r in data['records'])
    data['generated_at']=c.NOW
    data['schema_version']=2
    data['taxonomy']=c.TAXONOMY
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(data,ensure_ascii=False,separators=(',',':')));tmp.replace(path)
    print(json.dumps(data['relationship_scan'],indent=2),flush=True)

if __name__=='__main__':main()
