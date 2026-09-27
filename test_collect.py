import unittest
from collect import classify, record, matches, curated_geo_ids, collect_geo
from unittest.mock import patch

class ScientificClassification(unittest.TestCase):
    def test_reviewed_accessions_survive_empty_discovery(self):
        with patch('collect.CURATION', {'GSE123456': {}, 'E-MTAB-123': {}}), patch('collect.geo_search', return_value=([],0)), patch('collect.geo_fetch', return_value=[]) as fetch:
            self.assertEqual(curated_geo_ids(), ['200123456'])
            collect_geo(1)
            self.assertIn('200123456', fetch.call_args.args[0])

    def test_imaging_is_not_scrna(self):
        r=classify('Spatially resolved single-cell transcriptomic profiling of human lung adenocarcinoma using CosMx SMI')
        self.assertEqual(r['modalities'],['Spatial'])
        self.assertEqual(r['approaches'],['Imaging-based'])
    def test_nucleus_separate(self):
        r=classify('Pancreatic cancer with single-nucleus RNA-seq and Visium spatial transcriptomics')
        self.assertIn('snRNA-seq',r['modalities']);self.assertNotIn('scRNA-seq',r['modalities'])
    def test_no_donor_inference(self):
        r=record('test','Breast cancer with Visium and scRNA-seq','Matched computational analysis','https://example.org','GEO','2025-01-01')
        self.assertEqual(r['pairing'],'Both reported — unverified')
    def test_bulk_component_excluded(self):
        r=record('test','Cancer spatial atlas [bulk RNA-seq]','Visium and scRNA-seq in the study','https://example.org','GEO','2025-01-01')
        self.assertIsNone(r)
    def test_luad_is_included_in_nsclc_filter(self):
        self.assertIn('Non-small cell lung cancer',classify('Lung adenocarcinoma scRNA-seq')['subtypes'])
        self.assertNotIn('Non-small cell lung cancer',classify('Small cell lung cancer scRNA-seq')['subtypes'])

    def test_st_and_scrna_mentions_remain_unverified_study_leads(self):
        r=record('GSEtest','Lung adenocarcinoma [ST]','The study also generated scRNA-seq.','https://example.org','GEO','2021-11-30')
        self.assertEqual(set(r['modalities']),{'Spatial','scRNA-seq'})
        self.assertEqual(r['pairing'],'Both reported — unverified')

    def test_negated_subtype(self):
        self.assertFalse(matches('minimally invasive adenocarcinoma','invasive adenocarcinoma'))
        self.assertFalse(matches('non-small cell lung cancer','small cell lung cancer'))
        self.assertTrue(matches('small cell lung cancer','small cell lung cancer'))

if __name__=='__main__':unittest.main()
