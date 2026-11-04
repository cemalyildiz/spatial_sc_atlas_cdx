"""Evidence attached to repository samples; no donor-level data downloaded."""
import datetime as dt
import re
import urllib.request
import time
from urllib.parse import urlsplit,unquote
from collections import Counter

TISSUES={'Blood':r'peripheral blood|whole blood|\bPBMCs?\b|\bblood\b(?! vessel)','Intestine':r'\bgut\b|intestin|ileum|ileal|colon\b|colonic|rectum','Brain':r'brain|frontal cortex|prefrontal|cerebell|hippocamp','Lung':r'\blung\b|pulmonary','Kidney':r'kidney|renal cortex','Liver':r'\bliver\b','Heart':r'\bheart\b|myocard','Skin':r'\bskin\b|dermis','Synovium':r'synovi','Bone marrow':r'bone marrow','Breast':r'\bbreast\b','Lymph node':r'lymph node','Pancreas':r'pancrea','Spleen':r'spleen'}
DIAG={'Crohn’s disease':r'crohn','Ulcerative colitis':r'ulcerative colitis','Inflammatory bowel disease':r'\bIBD\b|inflammatory bowel disease','Alzheimer’s disease':r'alzheimer','Progressive supranuclear palsy (PSP)':r'progressive supranuclear palsy|\bPSP\b','Multiple sclerosis':r'multiple sclerosis','Rheumatoid arthritis':r'rheumatoid arthritis','Osteoarthritis':r'osteoarthritis','Lung adenocarcinoma':r'lung adenocarcinoma|\bLUAD\b','IgA nephropathy':r'IgA nephropathy'}
COHORT={'Pediatric':r'pediatric|paediatric','Adult':r'\badult','Treatment-naïve':r'treatment.na[iï]ve|untreated','Healthy control':r'healthy (?:donor|control|participant)|normal control','Control (unspecified)':r'\bcontrols?\b','Progressive disease':r'progressive disease','Postmortem':r'post.?mortem'}
FILE_TYPES={'Expression matrix':r'(?:matrix\.mtx|counts?\.h5|feature_bc_matrix|counts?_matrix)','Spatial coordinates':r'tissue_positions|spatial.coordinates|spot.coordinates','Cell boundaries':r'cell.boundar|cell.polygons|segmentation.mask','Histology image':r'tissue_(?:hires|lowres)_image|histology.*\.(?:tif|png|jpg)|h.?e.*\.svs','Annotations':r'cell.metadata|cell.annot|cell.type.*tsv','Sample map':r'sample.map|donor.map|sample.metadata'}
N='Not checked'
def match(text,terms):
    return sorted(k for k,v in terms.items() if re.search(v,text,re.I))
def download(a):
    url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+a+'&targ=gsm&form=text&view=brief'
    with urllib.request.urlopen(url,timeout=15) as res:raw=res.read(8_000_001)
    if len(raw)>8_000_000:raise ValueError('Response exceeds bounded sample metadata budget')
    samples={}
    for block in raw.decode().split('^SAMPLE = ')[1:]:
        lines=block.splitlines();sid=lines[0].strip()
        if not re.fullmatch(r'GSM\d+',sid):continue
        attrs=[];files=[]
        for line in lines:
            k,_,v=line.partition(' = ')
            if k=='!Sample_characteristics_ch1':
                key,sep,val=v.partition(':')
                if sep and re.fullmatch(r'tissue(?: type)?|organ|group|diagnosis|disease(?: state| status)?|condition|age group|treatment(?: status)?|control type|cell type',key.strip(),re.I):attrs.append([key.strip().lower(),val.strip()[:160]])
            elif k=='!Sample_source_name_ch1':attrs.append(['source',v[:160]])
            elif k=='!Sample_molecule_ch1':attrs.append(['molecule',v[:80]])
            elif k.startswith('!Sample_supplementary_file') and v.startswith(('ftp://ftp.ncbi.nlm.nih.gov/','https://ftp.ncbi.nlm.nih.gov/')):files.append(v.replace('ftp://','https://',1))
        samples[sid]={'attributes':attrs,'files':files}
    if not samples:raise ValueError('Missing samples')
    return samples

def profile(r):
    obs=[];files={};direct=r.get('gsm_facts',{})
    tissue=set();dx=set();clinical=set();mods=set();celltypes=set()
    for s in r.get('sample_metadata',[]):
        sid=s['id'];title=s['title'];a=direct.get(sid,{});attrs=a.get('attributes',[])
        celltypes.update(v for k,v in attrs if k=='cell type' and len(v)<100)
        t=' '.join([title]+[v for k,v in attrs if k in ('source','tissue','tissue type','organ')])
        c=' '.join([title]+[v for k,v in attrs if k not in ('tissue','organ','molecule','cell type')])
        ts=match(t,TISSUES);ds=match(c,DIAG);cs=match(c,COHORT)
        # Short forms are meaningful only in an explicitly identified series.
        if r['id']=='GSE350921' and any(k=='group' and re.search(r'\bAD\b',v) for k,v in attrs):ds+=['Alzheimer’s disease']
        if r['id']=='GSE350921' and any(k=='group' and 'control' in v.lower() for k,v in attrs):cs+=['Control (unspecified)']
        for k,v in attrs:
            if k in ('diagnosis','disease','disease state','disease status') and not match(v,DIAG) and 2<len(v)<80 and not re.search(r'normal|control|healthy|unknown|none|not applicable|negative|non.tumou?r|^no |^n/?a$',v,re.I):ds.append(v)
        if 'Healthy control' in cs and 'Control (unspecified)' in cs:cs.remove('Control (unspecified)')
        sampletypes=sorted({v for k,v in attrs if k=='cell type' and len(v)<100})
        if ts or ds or cs or sampletypes:
            obs.append(dict(id=sid,title=title,tissues=ts,diagnoses=sorted(set(ds)),clinical=sorted(set(cs)),cell_types=sampletypes,attributes=attrs))
            tissue.update(ts);dx.update(ds);clinical.update(cs)
        if any(k=='molecule' and v.lower()=='nuclear rna' for k,v in attrs) and re.search(r'single.nucleus RNA sequencing|snRNA.seq',r.get('title','')+' '+r.get('summary',''),re.I):mods.add('snRNA-seq')
        for url in a.get('files',[]):
            name=unquote(urlsplit(url).path.split('/')[-1])
            for k in match(name,FILE_TYPES) or ['Other listed files']:files.setdefault(k,[]).append({'name':name,'url':url,'sample':sid})
    if 'CELLxGENE' in r['sources']:
        for v in r.get('tissues',[]):tissue.update(match(v,TISSUES) or [v])
        for e in r.get('evidence',[]):
            if e.startswith('Repository disease annotations: '):dx.update(e.partition(': ')[2].split(', '))
        if r['access']=='Open processed data':
            for k in ['Expression matrix','Annotations']:files.setdefault(k,[]).append(dict(name='Repository collection & dataset downloads',url=r['url']))
    # Keep sample-deposited files separate from companion / controlled study resources.
    counts={};text=r.get('summary','')
    for k,pattern in {'participants':r'(?:Across|included|enrolled) ([\d,]+) participants','donors':r'from ([\d,]+) (?:unique )?(?:donors|patients)\b','cells':r'(?:contains|includes|comprises) ([\d,]+) (?:(?:peripheral|blood|immune|epithelial|tumou?r|B|T|single|human|quality.filtered) ){0,5}cells\b','nuclei':r'(?:contains|includes) ([\d,]+) nuclei','specimens':r'(?:sequencing of|from) ([\d,]+) [^.]{0,65}?specimens\b'}.items():
        m=re.search(pattern,text,re.I)
        if m:counts[k]={'n':int(m[1].replace(',','')),'evidence':m[0]}
    for k,pattern in {'donors':r'from ([\d,]+) (?:treatment.na.ve )?(?:patients|donors)\b','libraries':r'(?:across|as|by) (?:the )?([\d,]+) libraries'}.items():
        if k not in counts:
            m=re.search(pattern,r.get('overall_design','')+' '+text,re.I)
            if m:counts[k]={'n':int(m[1].replace(',','')),'evidence':m[0]}
    for url in r.get('series_files',[]):
        name=unquote(urlsplit(url).path.split('/')[-1])
        for k in match(name,FILE_TYPES) or ['Other listed files']:files.setdefault(k,[]).append({'name':name,'url':url})
    for rel in r.get('repository_relations',[]):
        m=re.fullmatch(r'SRA:\s*(https://(?:www\.)?ncbi\.nlm\.nih\.gov/\S+)',rel)
        if m:files.setdefault('Raw archive',[]).append(dict(name='SRA archive record · check access',url=m[1]))
    for l in r.get('links',[]):
        if 'controlled raw' in l['label'].lower() and l['url'].startswith('https://ega-archive.org/'):
            files.setdefault('Raw archive',[]).append(dict(name=l['label'],url=l['url'],status='Controlled access'))
    r.update(sample_cell_types=sorted(celltypes),sampled_tissues=sorted(tissue),sampled_diagnoses=sorted(dx),clinical_groups=sorted(clinical),sample_observations=obs,reported_counts=counts,file_inventory=files,file_categories=sorted(files),accession_modalities=sorted(set(r.get('accession_modalities',[]))|mods))

def refresh(records,old=(),limit=20):
    prev={r['id']:r for r in old}
    for r in records:
        for key in ('gsm_facts','gsm_checked','series_files'):
            if key not in r and key in prev.get(r['id'],{}):r[key]=prev[r['id']][key]
    priority=['GSE350921','GSE352630','GSE189487','GSE189357','GSE292469','GSE342824','GSE347958','GSE176078']
    todo=sorted((r for r in records if 'GEO' in r['sources'] and 0<(r.get('sample_count') or 0)<=200),key=lambda r:(r.get('gsm_checked',''),priority.index(r['id']) if r['id'] in priority else 999,r['id']))
    now=dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')
    for r in todo[:limit]:
        try:
            r['gsm_facts']=download(r['id'])
            time.sleep(.4)
            url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+r['id']+'&targ=self&form=text&view=brief'
            with urllib.request.urlopen(url,timeout=15) as res:body=res.read(2_000_000).decode()
            if '^SERIES = '+r['id'] not in body:raise ValueError('Series listing unavailable')
            r['series_files']=[v.replace('ftp://','https://',1) for v in re.findall(r'^!Series_supplementary_file = (\S+)',body,re.M) if v.startswith(('ftp://ftp.ncbi.nlm.nih.gov/','https://ftp.ncbi.nlm.nih.gov/'))]
            r['gsm_checked']=now
        except Exception as e:print('Sample metadata unavailable:',r['id'],str(e)[:120],flush=True)
        time.sleep(.4)
    for r in records:profile(r)
