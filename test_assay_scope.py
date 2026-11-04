import unittest
from assay_scope import scoped

def geo(title,**kwargs):
    return {'sources':['GEO'],'title':title,**kwargs}

class AssayScope(unittest.TestCase):
    def test_separate_xenium_not_inherited(self):
        r=geo('Human cancer [snRNA-seq]',summary='Matched Xenium measurements were deposited separately.',overall_design='Libraries were prepared with Chromium kits.')
        v=scoped(r)
        self.assertEqual(v['platforms'],['10x Chromium']);self.assertEqual(v['approaches'],['Sequencing-based'])
    def test_companion_and_reference_do_not_fill_unknown(self):
        for source in ['GEO','ArrayExpress']:
            v=scoped({'sources':[source],'title':'Human cancer atlas','summary':'Xenium and Visium are used in a companion study.','platforms':['Xenium'],'approaches':['Imaging-based']})
            self.assertEqual(v['platforms'],['Not specified']);self.assertEqual(v['approaches'],['Not specified'])
    def test_component_keeps_only_own_imaging_assay(self):
        v=scoped(geo('Lung cancer [Xenium in situ]',overall_design='Libraries were prepared using Chromium.'))
        self.assertEqual(v['platforms'],['Xenium']);self.assertEqual(v['approaches'],['Imaging-based'])
    def test_collection_ignores_summary_and_old_labels(self):
        v=scoped({'sources':['CELLxGENE'],'repository_assays':["10x 3' v3"],'summary':'Uses Xenium as reference','platforms':['Xenium']})
        self.assertEqual(v['platforms'],["10x 3' v3"]);self.assertEqual(v['approaches'],['Sequencing-based'])
    def test_related_records_and_sample_references_do_not_supply_assays(self):
        v=scoped(geo('A study',related=['Xenium'],sample_metadata=[{'title':'Visium external reference'},{'title':'CosMx tumor section'}]))
        self.assertEqual(v['platforms'],['CosMx'])
    def test_platform_does_not_imply_unspecified_readout(self):
        v=scoped(geo('Tumor [GeoMx]'))
        self.assertEqual(v['platforms'],['GeoMx']);self.assertEqual(v['approaches'],['Not specified'])
    def test_title_method(self):
        self.assertEqual(scoped(geo('Spatial transcriptomics of lung cancer using Visium HD'))['platforms'],['Visium HD'])

if __name__=='__main__':unittest.main()
