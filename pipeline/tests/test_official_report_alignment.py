"""Focused regression tests for scoped official metrics and report refresh."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scperturb.dragen.alignment import extract_scoped, extract_rank, display, prepare_alignment, SPECS, matches_official

class AlignmentTests(unittest.TestCase):
    def test_scopes_disambiguate_identical_labels(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'report.html'
            p.write_text('<div id="gex_grid"><script>rowData: [{"metric":"Input Reads","value":"100"}]</script></div><div id="all_grid"><script>rowData: [{"metric":"Input Reads","value":"130"}]</script></div>')
            out=extract_scoped(p)
            self.assertEqual(out['gex_grid','Input Reads']['value'],'100')
            self.assertEqual(out['all_grid','Input Reads']['value'],'130')
    def test_percentage_scale_and_count_format(self):
        self.assertEqual(display('96.00%'),'96.00%')
        self.assertEqual(display('77.9%'),'77.90%')
        self.assertEqual(display('168'),'168')
        self.assertEqual(display('11,112.5'),'11,112.5')
    def test_official_precision_and_exact_counts(self):
        self.assertTrue(matches_official(11112.5079,'11,112.5'))
        self.assertTrue(matches_official(45.692,'45.69%'))
        self.assertFalse(matches_official(45.70,'45.69%'))
        self.assertFalse(matches_official(168.1,'168'))
        self.assertTrue(matches_official(168,'168'))
    def test_duplicate_scoped_metric_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'report.html'; p.write_text('<div id="grid">rowData: [{"metric":"A","value":"1"},{"metric":"A","value":"2"}]')
            with self.assertRaises(ValueError): extract_scoped(p)
    def test_missing_or_mismatched_official_metric_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); (root/'03_Tables').mkdir(); (root/'03_Tables/project_summary.json').write_text('{"sample_id":"s"}')
            qc=root/'qc.json'; qc.write_text('{}')
            spec=SPECS[0]; _,scope,label,metric,*_=spec
            with patch('scperturb.dragen.alignment.expected_values',return_value=({metric:['','','','100']},{})),patch('scperturb.dragen.alignment.extract_scoped',return_value={(scope,label):{'value':'101'}}):
                with self.assertRaisesRegex(ValueError,'mismatch'): prepare_alignment(root,root,qc)
            with patch('scperturb.dragen.alignment.expected_values',return_value=({},{})),patch('scperturb.dragen.alignment.extract_scoped',return_value={}):
                with self.assertRaises(KeyError): prepare_alignment(root,root,qc)
    def test_official_rank_coordinates_preserved(self):
        traces=[{'label':'Cells','data':[{'x':1,'y':20},{'x':3,'y':5}]},{'label':'Background','data':[{'x':4,'y':4}]}]
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'r.html'; p.write_text('<div id="barcode_rank-subsection"><script>var rawTraces = '+json.dumps(traces)+';</script>')
            self.assertEqual(extract_rank(p),traces)
if __name__=='__main__': unittest.main()
