"""Review thresholds, sample isolation and report-only dispatch."""
import argparse
import csv
import gzip
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scperturb.validation import review
from scperturb import runner

ROOT=Path(__file__).resolve().parents[2]

class DeliveryReviewTests(unittest.TestCase):
    def test_threshold_boundaries_missing_and_nonfinite(self):
        self.assertEqual(review.evaluate(10000,">=",10000),"PASS")
        self.assertEqual(review.evaluate(9999,">=",10000),"NOT_PASS")
        self.assertEqual(review.evaluate(.60,"<=",.60),"PASS")
        self.assertEqual(review.evaluate(.60001,"<=",.60),"NOT_PASS")
        for value in (None,math.nan,math.inf):
            self.assertEqual(review.evaluate(value,">=",.25),"PENDING")
        self.assertEqual(review.evaluate(0,">=",.25,denominator=0),"PENDING")
        self.assertEqual(review.overall_result([{'status':'PASS'},{'status':'PENDING'}]),'PENDING')
        self.assertEqual(review.overall_result([{'status':'PENDING'},{'status':'NOT_PASS'}]),'NOT_PASS')

    def test_criteria_are_shared_with_standard(self):
        criteria=review.load_criteria(ROOT/review.DEFAULT_CRITERIA)
        self.assertEqual({r['id']:r['threshold'] for r in criteria['rules']}['QC01'],10000)
        self.assertEqual({r['id']:r['threshold'] for r in criteria['rules']}['QC02'],800)
        self.assertEqual(review.criteria_markdown(criteria),(ROOT/'pipeline/docs/SC_PERTURB_DELIVERY_REVIEW_CRITERIA.md').read_text())

    def test_invalid_percentage_threshold_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            c=json.loads((ROOT/review.DEFAULT_CRITERIA).read_text());c['rules'][2]['threshold']=80
            p=Path(d)/'criteria.json';p.write_text(json.dumps(c))
            with self.assertRaisesRegex(ValueError,'0..1'):review.load_criteria(p)

    def fixture(self,root,sample):
        import h5py
        def table(path,fields,rows):
            path.parent.mkdir(parents=True,exist_ok=True)
            with path.open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
        batch=root/'outputs'/sample; native=batch/'02_dragen_scrna_crispr/04_tool_output/dragen';native.mkdir(parents=True)
        delivery=batch/'05_report/04_tool_output'/('OUT_'+sample);tables=delivery/'03_Tables';tables.mkdir(parents=True)
        (root/'config').mkdir(); (root/'config/criteria.json').write_text((ROOT/review.DEFAULT_CRITERIA).read_text())
        (root/'config/analysis.json').write_text(json.dumps({'guide_reference':{'control_gene':'Non-Targeting Control'},'qc':{'min_target_cells':20},'reanalysis':{'guide_min_singlet_cells':5,'max_guide_fraction':.5}}))
        ref=root/'guides.csv';ref.write_text('gene,sgRNA_id,sequence\nTARGET_A,g1,AAAA\nTARGET_A,g2,CCCC\nNon-Targeting Control,n1,GGGG\nNon-Targeting Control,n2,TTTT\n')
        bc=[f'BC{i:03}' for i in range(120)]
        ids=['g1']*10+['g2']*10+['n1']*50+['n2']*50
        with (native/f'{sample}.scRNA.positive_cell_guide_assignments.csv').open('w') as f:
            w=csv.writer(f);w.writerow(['cell_barcode','num_features','feature_call','num_transcripts']);w.writerows((b,1,g,5) for b,g in zip(bc,ids))
        with (native/f'{sample}.scRNA.singlet_cell_guide_assignments.csv').open('w') as f:
            w=csv.writer(f);w.writerow(['cell_barcode','feature_call']);w.writerows(zip(bc,ids))
        with gzip.open(native/f'{sample}.scRNA.filtered.barcodes.tsv.gz','wt') as f:f.write('\n'.join(bc)+'\n')
        with gzip.open(native/f'{sample}.scRNA.features.tsv.gz','wt') as f:f.write('ENSG1\tG1\tGene Expression\nENSG2\tG2\tGene Expression\n')
        metrics={'Passing cells':120,'Total input reads':1200000,'Total input feature reads':96000,'Total barcoded reads':1080000,'Mapped reads':864000,'Fraction of reads in passing cells':.9,'Median genes per passing cell':1500,'Median molecules per passing cell':3000,'Total barcoded feature reads':86400,'Feature matching reads':69120,'Fraction of feature reads in passing cells':.9,'Fraction passing cells with CRISPR reads':1,'Filtered non-matching reads':10000}
        with (native/f'{sample}.scRNA_metrics.csv').open('w') as f:csv.writer(f).writerows(('TEST',sample,k,v) for k,v in metrics.items())
        with (native/f'{sample}.mapping_metrics.csv').open('w') as f:csv.writer(f).writerows([('MAPPING/ALIGNING SUMMARY','','Q30 bases',900,90),('MAPPING/ALIGNING SUMMARY','','Total input reads',1296000)])
        with (native/f'{sample}.scRNA.guide_metrics.csv').open('w') as f:csv.writer(f).writerows(('TEST',sample,k,v) for k,v in {'Total cells':120,'Total cells with one or more guides detected':120,'Total cells with two or more guides detected':0}.items())
        project={'sample_id':sample,'confidence_counts':{'PARTIAL_SUPPORT':1},'ntc_qc':{'ntc_singlet_cells':100,'active_ntc_guides_ge5_singlets':2,'ntc_guide_features':2},'dragen_native_crispr':True,'dragen_finished_normally':True,'dragen_license_success_in_source_log':True,'gex_reads':1200000,'crispr_feature_reads':96000,'sequencing_saturation':.2,'gsea_primary_tie_median':.8}
        (tables/'project_summary.json').write_text(json.dumps(project))
        table(tables/'dragen_official_metrics.tsv',['group','metric','display_value'],[{'group':'quality','metric':'% Q30','display_value':'90.00%'}])
        table(tables/'target_coverage_full.tsv',['target_gene','singlet_cells','active_guides_ge5','median_pairwise_spearman','median_directional_concordance'],[dict(target_gene='TARGET_A',singlet_cells=20,active_guides_ge5=2,median_pairwise_spearman=.1,median_directional_concordance=.6)])
        table(tables/'guide_coverage.tsv',['sgRNA_id','target_gene','singlet_cells'],[dict(sgRNA_id=g,target_gene='TARGET_A' if g.startswith('g') else 'Non-Targeting Control',singlet_cells=10 if g.startswith('g') else 50) for g in ['g1','g2','n1','n2']])
        table(tables/'target_confidence.tsv',['target_gene','confidence_class'],[dict(target_gene='TARGET_A',confidence_class='PARTIAL_SUPPORT')])
        table(tables/'screen_target_funnel.tsv',['count'],[{'count':1} for _ in range(5)])
        h5=batch/'03_cell_guide_qc/04_tool_output'/f'{sample}_gex_filtered.h5ad';h5.parent.mkdir(parents=True)
        with h5py.File(h5,'w') as f:
            obs=f.create_group('obs');obs.attrs['_index']='barcode';obs.create_dataset('barcode',data=bc,dtype=h5py.string_dtype())
            var=f.create_group('var');var.attrs['_index']='feature_id';var.create_dataset('feature_id',data=['ENSG1','ENSG2'],dtype=h5py.string_dtype())
            f.create_group('obsm').create_dataset('X_umap',shape=(120,2),dtype='f4')
        marker=batch/'02_dragen_scrna_crispr/05_metadata/stage_complete.json';marker.parent.mkdir(parents=True);marker.write_text(json.dumps({'status':'complete','sample_id':sample}))
        run={'project_root':str(root),'results_root':str(root/'outputs'),'guide_source':str(ref),'analysis_config':str(root/'config/analysis.json'),'delivery_review_config':str(root/'config/criteria.json')}
        return run,delivery

    def test_two_different_samples_and_read_requirements(self):
        for sample in ('alpha_new','beta_new'):
            with self.subTest(sample=sample),tempfile.TemporaryDirectory() as d:
                run,delivery=self.fixture(Path(d),sample)
                specs=[('gex','scope',f'M{i}','Total input reads',3,1,'input') for i in range(35)]
                scoped={('scope',f'M{i}'):{'value':'1200000'} for i in range(35)}
                with patch.object(review,'SPECS',specs),patch.object(review,'extract_scoped',return_value=scoped),patch.object(review,'verify_alignment',return_value=35):
                    result=review.collect_review(run,sample,delivery)
                review.apply_validation(result,{'result':'PASS','checks':[{'status':'PASS'}]*21})
                review.write_review(delivery,result)
                self.assertEqual(result['overall_result'],'PASS')
                self.assertEqual(result['observations']['gex_minimum_input_reads_at_current_passing_cells'],120*10000)
                self.assertEqual(result['observations']['crispr_minimum_input_reads_at_current_passing_cells'],120*800)
                self.assertTrue((delivery/f'04_Methods_Reproducibility/review/{sample}_DELIVERY_REVIEW.json').is_file())
                text=(delivery/f'04_Methods_Reproducibility/review/{sample}_DELIVERY_REVIEW.md').read_text()
                self.assertIn(sample,text)
                # A missed QC gate changes release status, not source values or analysis outputs.
                result['checks'][0]['status']='NOT_PASS';review.write_review(delivery,result)
                self.assertFalse(result['release_allowed']);self.assertEqual(result['overall_result'],'NOT_PASS')

    def test_cli_exit_tracks_technical_validation_not_qc_grade(self):
        with tempfile.TemporaryDirectory() as d:
            config=Path(d)/'run.json'
            config.write_text(json.dumps({'results_root':d}))
            result={'overall_result':'NOT_PASS','quantitative_qc_counts':{'NOT_PASS':1},'technical_counts':{'PASS':13}}
            for technical, expected in [('PASS',0),('FAIL',2)]:
                with patch.object(sys,'argv',['review','--config',str(config),'--sample-id','new_sample']), patch.object(review,'prepare_review',return_value=result), patch.object(review,'write_manifest'), patch.object(review,'validate',return_value={'result':technical}), patch.object(review,'finish_review',return_value=result):
                    self.assertEqual(review.main(),expected)

    def test_review_dispatches_no_analysis_stage(self):
        args=argparse.Namespace(config=Path('/project/config/run.yaml'))
        with patch.object(runner,'reports_ready',return_value=(Path('/project'),{},[{'sample_id':'different_sample'}])),patch.object(runner.subprocess,'run') as run,patch.object(runner,'verify'):
            runner.review(args)
        command=run.call_args.args[0]
        self.assertEqual(command[1:3],['-m','scperturb.validation.review'])
        self.assertIn('different_sample',command)

if __name__=='__main__':unittest.main()
