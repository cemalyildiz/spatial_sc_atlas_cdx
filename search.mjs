export const fields = ['sampled_tissues','sampled_diagnoses','clinical_groups','sample_cell_types','accession_modalities','file_categories','groups','diseases','subtypes','modalities','approaches','platforms','pairing','relationship_levels','sources','tissues','availability','access','validation'];
export function normalize(s) { return String(s ?? '').normalize('NFKD').replace(/[\u0300-\u036f]/g,'').replace(/[’‘]/g,"'").replace(/[–—]/g,'-').toLowerCase(); }
export function searchable(r) { return normalize([r.id,r.title,r.summary,...fields.flatMap(f=>r[f]||[]),...r.publications.flatMap(p=>[p.pmid,p.doi,p.title]),...(r.aliases||[]),...(r.related||[]),...(r.paired_sample_ids||[]),...(r.sample_metadata||[]).flatMap(s=>[s.id,s.title])].join(' ')); }
export function isBoth(r) { return r.modalities.includes('Spatial') && (r.modalities.includes('scRNA-seq') || r.modalities.includes('snRNA-seq')); }
export function priority(r) { return {'Matched donors / samples':5,'Same study — verified':4,'Same collection':3,'External reference':2,'Both reported — unverified':1}[r.pairing] || 0; }
export function filterRecords(records,state,skip) {
  const words=normalize(state.q||'').trim().split(/\s+/).filter(Boolean);
  return records.filter(r=>{
    if (words.length && !words.every(w=>(r._search||searchable(r)).includes(w))) return false;
    if (state.preset==='both'&&!isBoth(r)) return false;
    if (state.preset==='matched' && r.pairing!=='Matched donors / samples')return false;
    const local=state.companions?r.modalities:(r.accession_modalities||[]);
    if (state.preset==='spatial'&&!local.includes('Spatial')) return false;
    if (state.preset==='single'&&!local.some(m=>m==='scRNA-seq'||m==='snRNA-seq')) return false;
    if (state.from && (!r.released || r.released<state.from)) return false;
    if (state.to && (!r.released || r.released>state.to)) return false;
    const sampleFields=[['sampled_tissues','tissues'],['sampled_diagnoses','diagnoses'],['clinical_groups','clinical'],['sample_cell_types','cell_types']].filter(([f])=>f!==skip && state[f]?.length);
    if(sampleFields.length>1 && !(r.sample_observations||[]).some(s=>sampleFields.every(([f,key])=>state[f].some(v=>s[key]?.includes(v)))))return false;
    return fields.every(f=>(f==='platforms'&&state.unknownPlatforms&&r.platforms?.includes('Not specified')) || f===skip || !state[f]?.length || state[f].some(v=>(Array.isArray(r[f])?r[f]:[r[f]]).includes(v)));
  });
}
export function sortRecords(records,sort) {
  return [...records].sort((a,b)=>sort==='title'?a.title.localeCompare(b.title):sort==='oldest'? (a.released||'9999').localeCompare(b.released||'9999'):sort==='updated'?(b.checked||'').localeCompare(a.checked||''):sort==='newest'?(b.released||'').localeCompare(a.released||''):priority(b)-priority(a)||(b.released||'').localeCompare(a.released||''));
}
export function csvCell(value) {
  let text=Array.isArray(value)?value.join('; '):String(value??'');
  // Prevent spreadsheet formula evaluation when a repository title starts with =, +, -, @.
  if (/^[=+\-@\t\r]/.test(text)) text="'"+text;
  return '"'+text.replaceAll('"','""')+'"';
}
export function toCSV(records) {
  const columns=['id','title','sampled_tissues','sampled_diagnoses','clinical_groups','file_categories','diseases','subtypes','modalities','approaches','platforms','pairing','validation','relationship_levels','related','accession_modalities','paired_sample_count','sources','released','updated','checked','sample_count','cell_count','availability','access','url'];
  return '\uFEFF'+[columns.map(csvCell).join(','),...records.map(r=>columns.map(k=>csvCell(r[k])).join(','))].join('\r\n');
}
