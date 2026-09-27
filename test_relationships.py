import copy
import unittest
from relationships import components, build_relationships


def row(aid,mod,title='Human cancer atlas',design='',pmid='123',samples=()):
    return dict(id=aid,title=title,summary='',overall_design=design,accession_modalities=mod,modalities=mod[:],pairing='Not applicable',validation='Metadata-derived',evidence=[],publications=[{'pmid':pmid}],related=[],url='https://example.org/'+aid,sample_metadata=[{'title':s} for s in samples])

class Relationships(unittest.TestCase):
    def test_st_component_overrides_scrna_background(self):
        self.assertEqual(components('Spatiotemporal atlas of lung adenocarcinoma [ST]'),['Spatial'])
        self.assertEqual(components('Integrated spatial transcriptomics study [scRNA-seq]'),['scRNA-seq'])
    def test_shared_publication_and_sample_names_do_not_prove_matching(self):
        a=row('GSE1',['Spatial'],samples=['TD1 [ST]']);b=row('GSE2',['scRNA-seq'],samples=['TD1 scRNA-seq'])
        edges=build_relationships([a,b],{})
        self.assertEqual(edges[0]['status'],'Related study — unverified')
        self.assertEqual(edges[0]['shared_sample_labels'],['TD1'])
        self.assertIsNone(edges[0]['matched_sample_count'])
        self.assertEqual(a['modalities'],['Spatial'])
    def test_same_study_requires_title_design_and_publication(self):
        title='Spatiotemporal transcriptional atlas of human lung adenocarcinoma'
        design='Nine resected samples were collected for single-cell RNA sequencing and spatial transcriptomics.'
        a=row('GSE1',['Spatial'],title+' [ST]',design);b=row('GSE2',['scRNA-seq'],title+' [scRNA-seq]',design)
        edges=build_relationships([a,b],{})
        self.assertEqual(edges[0]['status'],'Same study — verified')
        self.assertEqual(a['accession_modalities'],['Spatial'])
        self.assertEqual(set(a['modalities']),{'Spatial','scRNA-seq'})
        self.assertFalse(edges[0]['reviewed'])
        self.assertEqual(build_relationships([a,b],{}),edges)
        b['overall_design']='A different cohort'
        self.assertEqual(build_relationships([a,b],{})[0]['status'],'Related study — unverified')
        self.assertEqual(a['modalities'],['Spatial'])
    def test_external_reference_requires_explicit_context(self):
        a=row('GSE1',['Spatial']);b=row('GSE2',['scRNA-seq'])
        a['summary']='Published scRNA-seq data GSE2 were used as an external reference.'
        self.assertEqual(build_relationships([a,b],{})[0]['status'],'External reference')
    def test_reviewed_subset_has_exact_count(self):
        a=row('GSE1',['Spatial']);b=row('GSE2',['scRNA-seq'])
        a.update(pairing='Matched donors / samples',evidence=['Six paired samples'],paired_sample_count=6,paired_sample_ids=['TD1','TD2','TD3','TD5','TD6','TD8'])
        edges=build_relationships([a,b],{'GSE1':{'related':['GSE2']}})
        self.assertEqual(edges[0]['matched_sample_count'],6)
        self.assertTrue(edges[0]['reviewed'])
    def test_doi_url_matches_bare_doi(self):
        a=row('GSE1',['Spatial']);b=row('GSE2',['scRNA-seq'])
        a['publications']=[{'doi':'https://doi.org/10.1234/ABC'}];b['publications']=[{'doi':'10.1234/abc'}]
        self.assertEqual(len(build_relationships([a,b],{})),1)

    def test_no_transitive_cohort_merge(self):
        a=row('GSE1',['Spatial'],pmid='1');b=row('GSE2',['scRNA-seq'],pmid='1');b['publications'].append({'pmid':'2'});c=row('GSE3',['Spatial'],pmid='2')
        edges=build_relationships([a,b,c],{})
        self.assertEqual({tuple(e['accessions']) for e in edges},{('GSE1','GSE2'),('GSE2','GSE3')})

if __name__=='__main__':unittest.main()
