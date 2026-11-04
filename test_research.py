import unittest,json
from pathlib import Path
from research import profile
from unittest.mock import patch
from collect import collect_geo

class ResearchEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows={x['id']:x for x in json.loads(Path('catalog.json').read_text())['records']}
    def test_blood_gut_not_eight_donors(self):
        r=self.rows['GSE352630']
        self.assertTrue({'Blood','Intestine'}<=set(r['sampled_tissues']))
        self.assertEqual(r['reported_counts']['participants']['n'],107)
        self.assertEqual(r['sample_count'],8)
        self.assertNotIn('donors',r['reported_counts'])
        self.assertTrue(any(x['id']=='GSM10159824' and 'Blood' in x['tissues'] for x in r['sample_observations']))
    def test_diagnosis_from_group_not_background(self):
        r=self.rows['GSE350921'];self.assertIn('Progressive supranuclear palsy (PSP)',r['sampled_diagnoses'])
        self.assertTrue(any('PSP' in v for s in r['sample_observations'] for k,v in s['attributes'] if k=='group'))
        self.assertFalse(any('corticobasal' in t.lower() for t in r['sampled_diagnoses']))
        self.assertIn('snRNA-seq',r['accession_modalities'])
    def test_counts_and_file_identity(self):
        r=self.rows['CXG:4f99649a-2021-47a3-9908-2b83c546cbc9']
        self.assertEqual(r['reported_counts']['cells']['n'],410646)
        self.assertEqual(r['reported_counts']['donors']['n'],37)
        spatial=self.rows['GSE189487'];self.assertTrue(all('GSM' in f['url'] for f in spatial['file_inventory']['Spatial coordinates']))
        self.assertNotIn('Spatial coordinates',self.rows['GSE189357']['file_inventory'])
    def test_missing_is_not_inferred(self):
        r=dict(id='TEST',sources=['GEO'],title='Brain cancer',summary='Other studies profile pediatric brain and blood with spatial images.',tissues=['Brain'])
        profile(r)
        self.assertEqual(r['sampled_tissues'],[]);self.assertEqual(r['file_categories'],[])
    def test_history_rotates_beyond_first_page(self):
        with patch('collect.TAXONOMY',[dict(disease='Example',terms=['example'],subtypes={})]),patch('collect.geo_search',return_value=(['1'],421)) as search,patch('collect.geo_fetch',return_value=[]):
            _,coverage=collect_geo(1)
            self.assertGreater(coverage['queries'][-2]['offset'],0)
            self.assertEqual(search.call_args.args[-1]%35,0)

if __name__=='__main__':unittest.main()
