"""Synthetic XF2 integration checks; no private recordings required."""
import contextlib
import io
import json
from pathlib import Path
import struct
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pyedflib

import decode_xf2


def fixture(path):
    blob = bytearray()
    for index in range(2):
        for kind, mapping, rate, values, dtype in (
            (160, 1, 4, [[0], [65535]], '<u2'),
            (161, 3, 4, [[-32768, 32767, -1, 0, 1, 2]] * 2, '>i2'),
        ):
            payload = np.array(values, dtype=dtype).tobytes()
            blob += b'\r' + struct.pack('<BIHH', kind, 1700000000, index * 500,
                                        6 + len(payload))
            blob += struct.pack('<HHHB', index, mapping, rate, 1)
            blob += payload + b'\0\0\n'
    path.write_bytes(blob)


class DecoderTests(unittest.TestCase):
    def run_cli(self, *args):
        with patch.object(sys, 'argv', ['decode_xf2.py', *map(str, args)]):
            with contextlib.redirect_stdout(io.StringIO()):
                decode_xf2.main()

    def test_file_both_and_folder_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'sample.XF2'
            fixture(source)
            self.run_cli(source, '--format', 'both')
            out = root / 'decoded'
            report = json.loads((out / 'conversion_report.json').read_text())[0]
            self.assertTrue(report['xdf_roundtrip_verified'])
            self.assertTrue(report['edf_roundtrip_verified'])
            with pyedflib.EdfReader(str(out / 'sample_EXG.edf')) as reader:
                np.testing.assert_array_equal(reader.readSignal(0), [0, 65535, 0, 65535])
            streams, _ = decode_xf2.decode(source)
            with np.load(out / 'sample_edf_timing.npz') as timing:
                for s in streams:
                    np.testing.assert_array_equal(timing[s['name'] + '_timestamps'], s['stamps'])
            self.run_cli(root, '--output-dir', root / 'other')
            self.assertTrue((root / 'other/sample.xdf').exists())
            self.assertFalse((root / 'other/sample_EXG.edf').exists())
            self.run_cli(source, '--format', 'edf', '--output-dir', root / 'edf_only')
            self.assertFalse((root / 'edf_only/sample.xdf').exists())

    def test_padding_and_scaling(self):
        with tempfile.TemporaryDirectory() as tmp:
            stream = dict(name='EXG', rate=4, labels=['Ch 01'],
                          values=np.array([[0], [32768], [65535]], dtype=np.int32),
                          stamps=1700000000 + np.arange(3) / 4)
            report = {}
            decode_xf2.write_edf(Path(tmp), 'short', [stream], report)
            self.assertEqual(report['edf']['streams']['EXG']['padding_samples'], 1)
            with pyedflib.EdfReader(str(Path(tmp) / 'short_EXG.edf')) as reader:
                np.testing.assert_array_equal(reader.readSignal(0), [0, 32768, 65535, 0])

    def test_invalid_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            for path in (Path(tmp), Path(tmp) / 'missing.xf2'):
                with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
                    self.run_cli(path)
                self.assertEqual(error.exception.code, 2)


if __name__ == '__main__':
    unittest.main()
