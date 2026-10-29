"""Reject incomplete or inconsistent benchmark evidence before publishing summaries."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

from validation.run_performance_audit import summarize


class AuditEvidence(unittest.TestCase):
    def fixture(self, directory, label='b1', *, seconds='0.2', signature='sig'):
        out = Path(directory)
        (out / f'{label}.manifest.json').write_text(json.dumps(dict(
            label=label, sizes=[2000], threads=1, reps=1, cases=['split'])))
        log = (f'AUDIT,{label},2000,split,1,0,0.3,0,2000,43\n'
               f'AUDIT,{label},2000,split,1,1,{seconds},0,2000,43\n'
               f'AUDIT_CHECK,{label},2000,split,1,{signature},.\n'
               f'AUDIT_COMPLETE label={label}\nCTOOLS_LAUNCH=abcd RC=0\n')
        (out / f'{label}.log').write_text(log)
        return out, log

    def reject(self, transform):
        with tempfile.TemporaryDirectory() as directory:
            out, log = self.fixture(directory)
            (out / 'b1.log').write_text(transform(log))
            with self.assertRaises(ValueError):
                summarize(out)
            self.assertFalse((out / 'summary.csv').exists())

    def test_complete(self):
        with tempfile.TemporaryDirectory() as directory:
            out, _ = self.fixture(directory)
            self.fixture(directory, 'c1')
            with contextlib.redirect_stdout(io.StringIO()):
                summarize(out)
            self.assertIn('b,2000,split,1,1,0.2', (out / 'summary.csv').read_text())
            self.assertEqual(len((out / 'raw.csv').read_text().splitlines()), 5)

    def test_missing_trial(self):
        self.reject(lambda log: log.replace('AUDIT,b1,2000,split,1,1,0.2,0,2000,43\n', ''))

    def test_missing_check(self):
        self.reject(lambda log: log.replace('AUDIT_CHECK,b1,2000,split,1,sig,.\n', ''))

    def test_duplicate_trial(self):
        self.reject(lambda log: log + 'AUDIT,b1,2000,split,1,1,0.2,0,2000,43\n')

    def test_duplicate_check(self):
        self.reject(lambda log: log + 'AUDIT_CHECK,b1,2000,split,1,sig,.\n')

    def test_missing_completion(self):
        self.reject(lambda log: log.replace('CTOOLS_LAUNCH=abcd RC=0', 'CTOOLS_LAUNCH=abcd RC=1'))

    def test_nonfinite_or_nonpositive(self):
        for value in ['nan', 'inf', '0', '-1']:
            with self.subTest(value=value):
                self.reject(lambda log: log.replace(',1,0.2,', f',1,{value},'))

    def test_different_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            out, _ = self.fixture(directory)
            self.fixture(directory, 'c1', signature='different')
            with self.assertRaisesRegex(ValueError, 'output mismatch'):
                summarize(out)
            self.assertFalse((out / 'summary.csv').exists())


if __name__ == '__main__':
    unittest.main()
