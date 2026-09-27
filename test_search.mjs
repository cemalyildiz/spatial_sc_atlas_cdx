import test from 'node:test';
import assert from 'node:assert/strict';
import {filterRecords,sortRecords,csvCell,isBoth} from './search.mjs';
const base={id:'a',title:'Study',summary:'',publications:[],released:'2026-09-01',checked:'2026-09-16',sources:['GEO'],subtypes:[],groups:['Cancer'],diseases:['Breast cancer'],modalities:['Spatial'],pairing:'Not applicable'};
const matched={...base,id:'b',released:'2021-09-01',subtypes:['Triple-negative (TNBC)'],modalities:['Spatial','scRNA-seq'],pairing:'Matched donors / samples'};
const nucleus={...base,id:'c',diseases:['Lung cancer'],modalities:['snRNA-seq']};
test('OR within disease, AND across facets',()=>{assert.deepEqual(filterRecords([base,matched,nucleus],{diseases:['Breast cancer','Lung cancer'],modalities:['snRNA-seq']}).map(r=>r.id),['c']);});
test('spatial-only is not paired, nucleus is not cell RNA',()=>{assert.equal(isBoth(base),false);assert.equal(isBoth(matched),true);assert.deepEqual(filterRecords([base,matched,nucleus],{preset:'both'}).map(r=>r.id),['b']);assert.equal(filterRecords([nucleus],{modalities:['scRNA-seq']}).length,0);});
test('unknown dates excluded when a date boundary is active',()=>{assert.equal(filterRecords([{...base,released:''},base],{from:'2026-01-01'}).length,1);});
test('source evidence outranks a newer unpaired record',()=>{assert.equal(sortRecords([base,matched],'priority')[0].id,'b');assert.equal(sortRecords([base,matched],'newest')[0].id,'a');});
test('subtype search and source accession search',()=>{assert.equal(filterRecords([base,matched],{q:'TNBC'}).length,1);assert.equal(filterRecords([base],{q:'nonexistent'}).length,0);});
test('CSV protects formula cells and preserves quotes',()=>{assert.equal(csvCell('=HYPERLINK("x")'),'"\'=HYPERLINK(""x"")"');assert.equal(csvCell('a,b'),'"a,b"');});

test('linked accession and relationship level are searchable and filterable',()=>{
  const linked={...base,related:['GSE189487'],relationship_levels:['Related study — unverified']};
  assert.equal(filterRecords([linked],{q:'GSE189487',relationship_levels:['Related study — unverified']}).length,1);
  assert.equal(filterRecords([linked],{relationship_levels:['Matched donors / samples']}).length,0);
});

test('reviewed LUAD pair is found with the original NSCLC and paired filters',async()=>{
  const {readFile}=await import('node:fs/promises');
  const catalogue=JSON.parse(await readFile(new URL('./catalog.json',import.meta.url),'utf8'));
  const result=filterRecords(catalogue.records,{q:'GSE189357',diseases:['Lung cancer'],subtypes:['Non-small cell lung cancer'],preset:'both',pairing:['Matched donors / samples']});
  assert.deepEqual(result.map(r=>r.id).sort(),['GSE189357','GSE189487']);
  assert.equal(result.find(r=>r.id==='GSE189357').sample_count,9);
  assert.equal(result.find(r=>r.id==='GSE189487').sample_count,6);
  assert.ok(result.every(r=>r.paired_sample_count===6));
});
