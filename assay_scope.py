"""Derive filter labels from record-local evidence, never from a linked study."""
import re

UNKNOWN = 'Not specified'
REFERENCE = re.compile(r'\breference\b|\bexternal\b|\bpreviously\b|\bseparately\b|\bcompanion\b|\bnot\b|\bwithout\b|\bcomparison\b|\bintegrat\w*\b|\breanaly\w*\b|\bGSE\d+', re.I)
RNA = re.compile(r'\b(?:s[cn]RNA[ -]?seq|single[ -](?:cell|nucleus) RNA[ -]?(?:seq|sequencing))',re.I)
SPATIAL = re.compile(r'spatial|\[ST\]|\b(?:Xenium|Visium|MERFISH|CosMx|GeoMx)\b',re.I)
CELL = {'10x Chromium','Smart-seq'}

def scoped(r, review=None):
    from collect import PLATFORMS
    evidence={}
    def note(platform,approach,kind,quote):
        key=(platform,approach)
        if key not in evidence:evidence[key]=dict(platform=platform,approach=approach,source=kind,text=quote)
    def readout(p,text):
        if p=='GeoMx' and not re.search(r'next.generation sequencing|\bNGS\b|Illumina',text,re.I):return None
        return PLATFORMS[p][1]
    def names(text):
        found={p for p,(pat,_) in PLATFORMS.items() if re.search(pat,text,re.I)}
        if '10x Chromium' in found and not re.search(r"Chromium|10x (?:Genomics )?[35][’'′]",text,re.I):
            found.remove('10x Chromium')
        if 'Visium HD' in found:found.discard('Visium')
        return found
    if 'CELLxGENE' in r['sources']:
        assays=r.get('repository_assays')
        if assays is None:
            # Older snapshots preserved the exact structured assay list in this evidence field.
            assays=next((s[len('Repository assays: '):].split(', ') for s in r.get('evidence',[]) if s.startswith('Repository assays: ')),[])
        for assay in assays:
            methods={readout(p,assay) for p in names(assay)}
            if not methods and re.search(r'RNA.seq|Smart.seq|Drop.seq|Seq.Well|10x [35]',assay,re.I):methods={'Sequencing-based'}
            for approach in methods or {None}:note(assay,approach,'Repository assay',assay)
    else:
        title=r['title']
        tags=re.findall(r'\[([^]]+)\]',title)
        suffix=re.search(r'\s[-–:]\s+([^:]+)$',title)
        if suffix:tags.append(suffix[1])
        scoped_title=' '.join(tags)
        exclusive_rna=bool(RNA.search(scoped_title)) and not SPATIAL.search(scoped_title)
        exclusive_spatial=bool(SPATIAL.search(scoped_title) or 'ST' in tags) and not RNA.search(scoped_title)
        def compatible(p):
            return not ((p in CELL and exclusive_spatial) or (p not in CELL and exclusive_rna))
        def use(text,kind):
            if REFERENCE.search(text):return
            for p in sorted(names(text)):
                if compatible(p):note(p,readout(p,text),kind,text)
            if RNA.search(text) and not exclusive_spatial:note(None,'Sequencing-based',kind,text)
        if not REFERENCE.search(title):
            for tag in tags:use(tag,'Component title')
            # A title naming one actual measurement method can identify its own assay.
            if len(names(title))==1 and not re.search(r'compar|benchmark|multi.platform|validation|reference|versus',title,re.I):
                p=next(iter(names(title)))
                if compatible(p) and not (RNA.search(title) and p not in CELL) and not (p in CELL and SPATIAL.search(title)):
                    if re.search(r'\b(?:using|with|by|on)\s+(?:the |10x (?:Genomics )?)?',title,re.I) or title.lower().startswith(p.lower()+' '):
                        use(title,'Record title')
        for sample in r.get('sample_metadata',[]):use(sample['title'],'Sample title')
        # Unqualified RNA-only titles support the approach, not a particular vendor.
        if RNA.search(title) and not SPATIAL.search(title) and not REFERENCE.search(title):
            note(None,'Sequencing-based','Record title',title)
        # Only explicit local production statements, with a matching component.
        # Skip mixed-method passages and background/reuse references.
        for field in ['overall_design','summary']:
            for line in re.split(r'(?<=[.!?])\s+',r.get(field,'') or ''):
                if REFERENCE.search(line) or len(names(line))!=1:continue
                p=next(iter(names(line)))
                if not ((p in CELL and exclusive_rna) or (p not in CELL and exclusive_spatial)):continue
                if re.search(r'\b(?:libraries|samples|sections|specimens|nuclei|cells)\b.{0,65}\b(?:prepared|processed|profiled|sequenced|assayed|analy[sz]ed)\b.{0,25}\b(?:using|with|on)\b',line,re.I):
                    use(line,'Local method statement')
        if review:
            for p in review['platforms']:note(p,readout(p,review['basis']),'Reviewed accession: '+review['url'],review['basis'])
    platforms=sorted({p for p,_ in evidence if p})
    approaches=sorted({a for _,a in evidence if a})
    return dict(platforms=platforms or [UNKNOWN],approaches=approaches or [UNKNOWN],assay_evidence=list(evidence.values()))

def audit(records,curation):
    for r in records:r.update(scoped(r,curation.get(r['id'],{}).get('accession_assay_review')))
