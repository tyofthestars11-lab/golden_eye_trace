"""Filterbank reader/writer tests, including markers seen in real files."""

import os
import struct
import tempfile
import unittest

import numpy as np

from seti_search.filterbank import Filterbank, parse_header, write_filterbank


def _string(value: str) -> bytes:
    raw = value.encode()
    return struct.pack("<i", len(raw)) + raw


class ReaderTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def path(self, name: str = "test.fil") -> str:
        return os.path.join(self.tmp.name, name)

    def test_float32_roundtrip(self):
        rng = np.random.default_rng(0)
        data = rng.exponential(size=(32, 256)).astype(np.float32)
        path = self.path()
        write_filterbank(path, data, fch1=1420.5, foff=-2.86e-6, tsamp=0.5, source="TEST")

        fb = Filterbank.open(path)
        self.assertEqual(fb.nchans, 256)
        self.assertEqual(fb.nints, 32)
        self.assertEqual(fb.nbits, 32)
        self.assertEqual(fb.source, "TEST")
        self.assertAlmostEqual(fb.tsamp, 0.5)
        np.testing.assert_allclose(fb.read(), data, rtol=1e-6)

        freqs = fb.frequencies()
        self.assertAlmostEqual(freqs[0], 1420.5)
        self.assertAlmostEqual(freqs[-1], 1420.5 - 2.86e-6 * 255, places=9)
        # foff is negative, so the stored channel order descends in frequency
        self.assertGreater(freqs[0], freqs[-1])

    def test_uint8_roundtrip(self):
        data = (np.arange(10 * 64) % 251).reshape(10, 64).astype(np.uint8)
        path = self.path("eight.fil")
        write_filterbank(path, data, fch1=1400.0, foff=0.001, tsamp=2.0, nbits=8)
        fb = Filterbank.open(path)
        self.assertEqual(fb.nbits, 8)
        np.testing.assert_array_equal(fb.read(), data.astype(np.float32))

    def test_subrange_read(self):
        data = np.arange(6 * 8, dtype=np.float32).reshape(6, 8)
        path = self.path("sub.fil")
        write_filterbank(path, data, fch1=100.0, foff=0.5, tsamp=1.0)
        fb = Filterbank.open(path)
        np.testing.assert_array_equal(fb.read(2, 5, 3, 7), data[2:5, 3:7])

    def test_raw_markers(self):
        """Some writers emit HEADER_START / HEADER_END as raw bytes."""
        nchans, nints = 4, 3
        data = np.arange(nints * nchans, dtype=np.float32).reshape(nints, nchans)

        header = b"HEADER_START" + struct.pack("<i", 3093750)
        header += _string("RAW")
        header += struct.pack("<i", 6)  # telescope_id
        header += struct.pack("<i", 0)  # machine_id
        header += struct.pack("<i", 1)  # data_type
        header += struct.pack("<i", 0)  # barycentric
        header += struct.pack("<i", 0)  # pulsarcentric
        for value in (0.0, 0.0, 0.0, 0.0):  # az, za, raj, dej
            header += struct.pack("<d", value)
        header += struct.pack("<d", 0.0)  # tstart
        header += struct.pack("<d", 1.0)  # tsamp
        header += struct.pack("<i", 32)  # nbits
        header += struct.pack("<d", 100.0)  # fch1
        header += struct.pack("<d", 0.25)  # foff
        header += struct.pack("<i", nchans)
        header += struct.pack("<i", 1)  # nifs
        header += b"HEADER_END"

        path = self.path("raw.fil")
        with open(path, "wb") as handle:
            handle.write(header)
            handle.write(data.astype("<f4").tobytes())

        fb = Filterbank.open(path)
        self.assertEqual(fb.source, "RAW")
        self.assertEqual(fb.nchans, nchans)
        self.assertEqual(fb.nints, nints)
        self.assertAlmostEqual(fb.foff, 0.25)
        np.testing.assert_array_equal(fb.read(), data)

    def test_extra_keys_are_typed(self):
        """nbeams/ibeam appear in multibeam files before source_name."""
        header = _string("HEADER_START") + struct.pack("<i", 3093750)
        header += _string("nbeams") + struct.pack("<i", 13)
        header += _string("ibeam") + struct.pack("<i", 7)
        header += _string("HEADER_END")
        parsed, _ = parse_header(header)
        self.assertEqual(parsed["nbeams"], 13)
        self.assertEqual(parsed["ibeam"], 7)

    def test_rejects_unknown_key(self):
        header = _string("HEADER_START") + struct.pack("<i", 3093750)
        header += _string("mystery_key") + struct.pack("<i", 1)
        header += _string("HEADER_END")
        with self.assertRaises(ValueError) as caught:
            parse_header(header)
        self.assertIn("mystery_key", str(caught.exception))

    def test_rejects_non_filterbank(self):
        path = self.path("nope.fil")
        with open(path, "wb") as handle:
            handle.write(b"<html><body>this is a web page</body></html>")
        with self.assertRaises(ValueError):
            Filterbank.open(path)

    def test_partial_trailing_record_is_ignored(self):
        data = np.ones((4, 16), dtype=np.float32)
        path = self.path("partial.fil")
        write_filterbank(path, data, fch1=1.0, foff=0.1, tsamp=1.0)
        with open(path, "ab") as handle:
            handle.write(b"\x00" * 12)  # less than one full integration
        fb = Filterbank.open(path)
        self.assertEqual(fb.nints, 4)


if __name__ == "__main__":
    unittest.main()
