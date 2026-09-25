export const fields = ['groups','diseases','subtypes','modalities','approaches','platforms','pairing','sources','tissues','availability','access','validation'];
export function normalize(s) { return String(s ?? '').normalize('NFKD').replace(/[\u0300-\u036f]/g,'').replace(/[’‘]/g,"'").replace(/[–—]/g,'-').toLowerCase(); }
export function searchable(r) { return normalize([r.id,r.title,r.summary,...fields.flatMap(f=>r[f]||[]),...r.publications.flatMap(p=>[p.pmid,p.doi,p.title]),...(r.aliases||[])].join(' ')); }
export function isBoth(r) { return r.modalities.includes('Spatial') && (r.modalities.includes('scRNA-seq') || r.modalities.includes('snRNA-seq')); }
export function priority(r) { return {'Matched donors / samples':5,'Same study — verified':4,'Same collection':3,'External reference':2,'Both reported — unverified':1}[r.pairing] || 0; }
export function filterRecords(records,state,skip) {
  const words=normalize(state.q||'').trim().split(/\s+/).filter(Boolean);
  return records.filter(r=>{
    if (words.length && !words.every(w=>(r._search||searchable(r)).includes(w))) return false;
    if (state.preset==='both'&&!isBoth(r)) return false;
    if (state.preset==='spatial'&&!r.modalities.includes('Spatial')) return false;
    if (state.preset==='single'&&!r.modalities.some(m=>m==='scRNA-seq'||m==='snRNA-seq')) return false;
    if (state.from && (!r.released || r.released<state.from)) return false;
    if (state.to && (!r.released || r.released>state.to)) return false;
    return fields.every(f=>f===skip || !state[f]?.length || state[f].some(v=>(Array.isArray(r[f])?r[f]:[r[f]]).includes(v)));
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
  const columns=['id','title','diseases','subtypes','modalities','approaches','platforms','pairing','validation','sources','released','updated','checked','sample_count','cell_count','availability','access','url'];
  return '\uFEFF'+[columns.map(csvCell).join(','),...records.map(r=>columns.map(k=>csvCell(r[k])).join(','))].join('\r\n');
}
