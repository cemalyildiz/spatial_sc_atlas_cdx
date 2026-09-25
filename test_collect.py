import unittest
from collect import classify, record, matches

class ScientificClassification(unittest.TestCase):
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
    def test_negated_subtype(self):
        self.assertFalse(matches('non-small cell lung cancer','small cell lung cancer'))
        self.assertTrue(matches('small cell lung cancer','small cell lung cancer'))

if __name__=='__main__':unittest.main()
