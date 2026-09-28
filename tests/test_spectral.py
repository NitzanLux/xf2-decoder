import unittest
import numpy as np
from spectral_analysis import spectra


class SpectralTests(unittest.TestCase):
    def test_identical_signals_and_known_frequency(self):
        fs = 500
        t = np.arange(5000)/fs
        x = (np.sqrt(2)*np.sin(2*np.pi*40*t))[:, None]
        result, meta = spectra(x, x, fs)
        self.assertEqual(result['frequency_hz'][result['lsl_psd'][:, 0].argmax()], 40)
        self.assertAlmostEqual(np.sum(result['lsl_psd'][:, 0])*0.5, 1, places=5)
        np.testing.assert_allclose(result['coherence'], 1, atol=1e-12)
        np.testing.assert_array_equal(result['residual_psd'], 0)
        self.assertEqual(meta['averaged_segments'], 9)

    def test_independent_noise_and_zero_power(self):
        rng = np.random.default_rng(42)
        result, _ = spectra(rng.normal(size=(5000, 4)), rng.normal(size=(5000, 4)), 500)
        self.assertLess(np.mean(result['coherence']), .2)
        empty, _ = spectra(np.zeros((5000, 1)), np.zeros((5000, 1)), 500)
        self.assertTrue(np.isnan(empty['coherence']).all())

    def test_too_short(self):
        with self.assertRaises(ValueError):
            spectra(np.ones((1000, 1)), np.ones((1000, 1)), 500)


if __name__ == '__main__':
    unittest.main()
