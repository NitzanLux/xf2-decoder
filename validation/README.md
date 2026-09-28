# Validation: LSL versus decoded SD

A reproducible comparison of LSL recordings against the corresponding SD recordings decoded by [`decode_xf2.py`](../decode_xf2.py), over **LSL sample-index time 100–110 seconds**. It includes small original-value recording excerpts, analysis code, per-channel metrics, and generated figures. No full recordings or external data paths are required.

Run all commands from the repository root, after following the install steps in the [main README](../README.md#install).

## Layout

| Path | Contents |
| --- | --- |
| [`compare_exg.py`](compare_exg.py) | 16-channel EXG alignment, metrics, and plots |
| [`compare_imu.py`](compare_imu.py) | 6-axis IMU alignment, metrics, and plots |
| [`spectral_analysis.py`](spectral_analysis.py) | Welch PSD, residual PSD, and coherence for both |
| [`data/`](data/) | Cropped raw excerpts and source provenance |
| [`results/exg/`](results/exg/) | EXG outputs |
| [`results/imu/`](results/imu/) | IMU outputs |

Both results folders use the same filenames:

| File | Contents |
| --- | --- |
| `aligned_all_channels_10s.png` / `.svg` | Aligned standardized traces with per-channel scores |
| `residual_all_channels_10s.png` / `.svg` | Standardized LSL minus SD, common vertical scale |
| `channel_metrics.csv` | Per-channel Pearson r, standardized RMSE (EXG also: r², fitted gain, fit RMSE) |
| `summary.json` | Alignment offset and aggregate metrics |
| `comparison_arrays.npz` | Aligned standardized samples and per-sample `residual_z` |
| `psd_10s`, `residual_psd_10s`, `coherence_10s` (`.png` / `.svg`) | Spectral plots |
| `spectral_arrays.npz`, `spectral_metadata.json` | Spectral values and settings |

`results/exg/` also contains `convolution_crosscorrelation.png`, and its `comparison_arrays.npz` includes the full convolution and cross-correlation arrays.

## EXG comparison

```sh
.venv/bin/python validation/compare_exg.py
```

![All aligned channels](results/exg/aligned_all_channels_10s.png)

Open [the SVG](results/exg/aligned_all_channels_10s.svg) for lossless zoom, or [the CSV](results/exg/channel_metrics.csv) for quantitative results. Each channel title reports Pearson r and standardized RMSE for the actual displayed 10 seconds. The separate [residual plot (SVG)](results/exg/residual_all_channels_10s.svg) ([PNG](results/exg/residual_all_channels_10s.png)) shows standardized LSL minus standardized SD, with a shared vertical scale across channels. Zero means the standardized traces agree at that sample; positive values mean LSL is higher. Residuals are in z-score units, not µV or raw ADC counts.

The SD curve uses alpha 0.45 and LSL uses 0.95, revealing both overlapping traces. Alpha is **opacity**: lower values give more transparency. To increase SD opacity:

```sh
.venv/bin/python validation/compare_exg.py --alpha 0.9
```

### Included EXG data

`data/raw_sample.npz` contains only 18 seconds of LSL data (96–114 s, 500 Hz, µV) and 20 seconds of SD data (150–170 s, 4,000 Hz, raw counts). Padding supports filtering and alignment of the central 10-second example. The arrays are `lsl` and `sd`, both samples × 16 channels. Source filenames, SHA-256 hashes, sample ranges, units, and rates are recorded in `data/provenance.json`. Times refer to sample indices divided by nominal rates, not synchronized absolute XDF timestamps.

### EXG processing and interpretation

1. Anti-alias FIR decimation reduces SD from 4,000 Hz to 500 Hz using zero-phase decimation by 8.
2. Both signals receive a fourth-order 20 Hz Butterworth high-pass applied forward and backward. Padding reduces boundary effects.
3. Normalized cross-correlation searches the excerpt for one shared offset across channels. It uses an interior LSL segment, excluding two seconds at each edge.
4. The aligned 100–110 s interval is extracted. Each channel and recording is standardized independently to zero mean and unit population standard deviation over those 10 seconds.
5. Pearson r, r², standardized RMSE, and a linear fit from LSL µV to SD counts are calculated per channel. The linear fit's RMSE is in counts. Fitted gains are empirical comparisons, not physical calibration.

Standardization supports waveform comparison but removes absolute amplitude differences; consult the fitted gains and residuals in the CSV for that aspect. Alignment and scores use overlapping data, so these are descriptive scores, not independent validation. This excerpt lies within one stable-offset segment. The full recording has several timing steps and needs local alignment; one offset must not be assumed for its entire duration. Cropped filtering may produce tiny numerical differences from filtering the complete files.

### What causes the residual error?

The [residual plot](results/exg/residual_all_channels_10s.svg) shows **standardized LSL minus standardized SD** after filtering and alignment. Zero means the standardized traces agree at that sample. These differences are in z-score units, not µV or raw ADC counts; they measure waveform agreement after each trace's mean and amplitude scale have been removed.

**Downsampling can contribute, but this comparison does not isolate or quantify its contribution.** Reducing SD from 4,000 Hz to 500 Hz applies an anti-alias FIR filter before decimation, attenuating high-frequency content. If the LSL acquisition path uses a different filter response, the resulting waveforms can differ. Downsampling does not necessarily cause a mismatch when both paths have compatible bandwidth and timing. The LSL acquisition filter response has not been established here.

Other possible contributors include a remaining fractional-sample timing offset, clock drift within the window, differences in acquisition filtering, and noise. The current alignment uses one shared integer-sample shift; at 500 Hz, one sample is 2 ms. Even a smaller offset can produce visible residuals around sharp peaks. Independent standardization also affects the displayed differences, so this plot does not measure absolute amplitude accuracy.

To assess the downsampling contribution, first establish the acquisition filter responses, compare the signals over a matched bandwidth, and check how fractional-sample alignment and local clock correction change the residual. Those checks have not been performed in this example. **The residual alone does not establish a decoder error or show that downsampling is its main cause.**

### Convolution and cross-correlation

[`results/exg/convolution_crosscorrelation.png`](results/exg/convolution_crosscorrelation.png) shows curves for each channel, zoomed to ±100 ms. The curve labeled **Time-reversed overlap (convolution)** describes multiplication and summation with one signal reversed in time. **Time-shift similarity (cross-correlation)** measures agreement at different shifts. Convolution and cross-correlation operate on standardized 10-second signals and are divided by N. Nonzero-lag correlation uses fixed-N normalization, not overlap-specific Pearson normalization. The convolution plot centers its axis for display; the saved array includes the true output time axis. Convolution is not itself a similarity score.

## IMU comparison

```sh
.venv/bin/python validation/compare_imu.py
```

Open the [aligned IMU traces (SVG)](results/imu/aligned_all_channels_10s.svg) ([PNG](results/imu/aligned_all_channels_10s.png)) and the separate [IMU residual plot (SVG)](results/imu/residual_all_channels_10s.svg) ([PNG](results/imu/residual_all_channels_10s.png)). Both show Acc X/Y/Z and Gyro X/Y/Z over LSL IMU sample-index time 100–110 seconds. [Per-axis metrics](results/imu/channel_metrics.csv) report Pearson r and standardized RMSE.

SD IMU is decimated from 1,000 Hz to 500 Hz with a zero-phase anti-alias FIR filter. **No 20 Hz high-pass filter is applied to IMU**, preserving slow motion and orientation changes. A shared integer-sample offset across all six axes is estimated independently of EXG by averaging normalized cross-correlations over the central LSL excerpt. The selected offset is 55.102 seconds in nominal sample-index time; it is not an absolute timestamp correction and should not be applied to EXG or other recordings.

Each source and axis is standardized separately over the displayed window. LSL header units are g for acceleration and deg/s for angular velocity; decoded SD values are raw counts without established physical calibration. Consequently, these plots compare waveform shape, not physical amplitude accuracy. Residuals are standardized LSL minus standardized SD, with a common vertical scale across axes. Filtering, remaining timing differences, and noise can contribute to the residual, as discussed above for EXG.

For this excerpt, Pearson r ranges from 0.99930 to 0.99990 across the six axes. Alignment and evaluation use overlapping data, so this is a descriptive comparison, not independent validation. The bundled `data/imu_raw_sample.npz` retains original values for LSL 96–114 seconds and SD 145–180 seconds from the same source XDF files used for EXG. `data/imu_provenance.json` records hashes, stream names, sample ranges, labels, and source units. The example runs without the full recordings. `results/imu/comparison_arrays.npz` contains aligned standardized samples and residuals; `results/imu/summary.json` records the alignment and aggregate metrics.

## Spectral analysis: EXG and IMU

Run spectral analysis on the saved aligned, standardized 100–110 second segments:

```sh
# Refresh the aligned arrays if you change comparison processing:
.venv/bin/python validation/compare_exg.py
.venv/bin/python validation/compare_imu.py
.venv/bin/python validation/spectral_analysis.py
```

Use `--signal exg` or `--signal imu` to analyze only one signal type. The default is both.

| Analysis | EXG | IMU |
| --- | --- | --- |
| LSL and SD power spectral density | [PSD](results/exg/psd_10s.svg) | [PSD](results/imu/psd_10s.svg) |
| Residual power spectral density | [Residual PSD](results/exg/residual_psd_10s.svg) | [Residual PSD](results/imu/residual_psd_10s.svg) |
| Magnitude-squared coherence | [Coherence](results/exg/coherence_10s.svg) | [Coherence](results/imu/coherence_10s.svg) |

PNG versions are also generated. Power spectral density (PSD) shows how signal power is distributed across frequencies. Similar LSL and SD curves indicate similar spectral shapes. The separate residual PSD is computed from **LSL z-score minus SD z-score**, not by subtracting their power spectra; it shows which frequencies contribute to the remaining disagreement.

Coherence measures the consistency of the linear relationship between the recordings at each frequency, on a scale from 0 to 1. A value near 1 does not establish equal amplitudes, zero delay, or correct calibration. Interpret it alongside PSD: estimates where signal power is very small are less informative. Zero-power bins have undefined coherence and are saved as NaN.

The implementation uses [SciPy Welch spectral estimation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html), with two-second Hann windows, 50% overlap, mean removal within each window, and nine averaged segments. At 500 Hz this gives 0.5 Hz frequency-bin spacing; the window limits the ability to resolve nearby components. The short, overlapping segments provide a descriptive estimate, not a statistical significance test.

All PSD units are **z-score²/Hz**, because each displayed trace was standardized independently. These are comparisons of relative spectral content, not µV²/Hz, g²/Hz, or physical noise floors. EXG already has the 20 Hz high-pass filter applied; its shaded region below 20 Hz is attenuated by preprocessing. IMU has no high-pass filter and uses a logarithmic frequency axis to show slower motion. Both cover positive frequencies through the 250 Hz Nyquist limit; DC is retained in the saved arrays but omitted from plots. Anti-alias filtering affects the upper frequencies, and frequencies above 250 Hz in the original SD recordings are not analyzed here.

`results/exg/spectral_arrays.npz` and `results/imu/spectral_arrays.npz` store frequency bins, both PSDs, residual PSD, and coherence. The corresponding `spectral_metadata.json` files record the input array file, window settings, interval, units, and preprocessing. Spectral estimates reuse the existing alignment; they do not independently validate the decoder or identify the cause of a discrepancy.
