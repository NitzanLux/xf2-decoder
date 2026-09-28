# XF2 decoder with XDF and EDF export

[`decode_xf2.py`](decode_xf2.py) is the main tool: it decodes the observed XtrodES XF2 layout into XDF or EDF+, preserving raw ADC counts and original packet timing. The included LSL-versus-SD comparison provides a small, reproducible validation example.

## Install

Python 3.9 or later:

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Decode XF2 recordings

Pass a folder containing `.xf2` files:

```sh
.venv/bin/python decode_xf2.py "/path/to/recordings"
```

The decoder processes every `.xf2` file directly in that folder and writes into its `decoded/` subfolder:

- One `.xdf` per source file, containing available EXG/IMU signals and corresponding packet-timing streams.
- Per-stream PNG trace plots in raw counts.
- `conversion_report.json` with source hashes, packet continuity, clock-fit statistics, and XDF round-trip verification.
- `README.txt` describing the export.

Each exported XDF is reloaded with `pyxdf`; sample values and timestamps are checked for exact equality against the decoded arrays. The decoder checks observed frame boundaries, payload sizes, and packet continuity. It rejects missing packets rather than silently filling gaps. It supports the observed format and constant stream configuration; it is not a general specification of all XF2 variants. CRC bytes are present but their algorithm is unverified. Physical calibration is not established, and device timestamps are not synchronized to LSL host time.

Pass a single file instead of a folder, and select EDF or both formats:

```sh
.venv/bin/python decode_xf2.py "/path/to/recording.xf2" --format edf
.venv/bin/python decode_xf2.py "/path/to/recordings" --format both --output-dir "/path/to/exports"
```

XDF remains the default. A single-file input writes to `decoded/` beside that file unless `--output-dir` is supplied. Folder input is nonrecursive and accepts case-insensitive `.xf2` extensions. Invalid paths and folders without XF2 files produce an error.

EDF export creates `<stem>_EXG.edf` and/or `<stem>_IMU.edf`, one EDF+ file per available signal stream, so separate start times and lengths are retained. EDF uses each stream's **nominal integer sample rate**, not its fitted clock. Physical values are uncalibrated raw counts: unsigned EXG counts are stored with a digital offset of -32768 and mapped back through the EDF physical range; IMU counts remain signed. Every exported channel is reloaded to verify digital samples, physical count scaling, and sample rate.

EDF requires complete records: the last one-second record is padded with raw zero counts and annotated. `conversion_report.json` records original sample counts and padding. `<stem>_edf_timing.npz` preserves exact fitted timestamps (`EXG_timestamps`, `IMU_timestamps`) and packet timing arrays (`EXG_packet_timing_timestamps`, `EXG_packet_timing_values`, and IMU equivalents). Packet value columns are packet index, first sample index, and sample count. Keep these sidecars and the report with the EDF files. EDF start dates interpret device Unix seconds as UTC; this does not imply clock synchronization. EDF export requires integer nominal rates. This tool reads XF2 input; it does not reconstruct XF2 from EDF.

Existing files with matching output names are overwritten. Run with normal Python, without `-O`, so the XF2 and XDF assertion checks remain enabled.

No full XF2 recordings are bundled. Supply your own folder for decoding; the validation example below runs entirely from the included cropped samples.

## LSL-versus-decoded-SD validation example

A reproducible comparison of 16 LSL channels and their corresponding SD channels over **LSL sample-index time 100–110 seconds**. Includes small original-value recording excerpts, analysis code, per-channel metrics, and generated figures. No full recordings or external data paths are required.

![All aligned channels](results/aligned_all_channels_10s.png)

Open [the SVG](results/aligned_all_channels_10s.svg) for lossless zoom, or [the CSV](results/channel_metrics.csv) for quantitative results. Each channel title reports Pearson r and standardized RMSE for the actual displayed 10 seconds. The separate [residual plot (SVG)](results/residual_all_channels_10s.svg) ([PNG](results/residual_all_channels_10s.png)) shows standardized LSL minus standardized SD, with a shared vertical scale across channels. Zero means the standardized traces agree at that sample; positive values mean LSL is higher. Residuals are in z-score units, not µV or raw ADC counts.

## Run the comparison

After installing the dependencies:

```sh
.venv/bin/python compare.py
```

The SD curve uses alpha 0.45 and LSL uses 0.95, revealing both overlapping traces. Alpha is **opacity**: lower values give more transparency. To increase SD opacity:

```sh
.venv/bin/python compare.py --alpha 0.9
```

## IMU traces comparison

Run the six-axis accelerometer/gyroscope comparison from the bundled IMU excerpts:

```sh
.venv/bin/python compare_imu.py
```

Open the [aligned IMU traces (SVG)](results/imu_aligned_all_channels_10s.svg) ([PNG](results/imu_aligned_all_channels_10s.png)) and the separate [IMU residual plot (SVG)](results/imu_residual_all_channels_10s.svg) ([PNG](results/imu_residual_all_channels_10s.png)). Both show Acc X/Y/Z and Gyro X/Y/Z over LSL IMU sample-index time 100–110 seconds. [Per-axis metrics](results/imu_channel_metrics.csv) report Pearson r and standardized RMSE.

SD IMU is decimated from 1,000 Hz to 500 Hz with a zero-phase anti-alias FIR filter. **No 20 Hz high-pass filter is applied to IMU**, preserving slow motion and orientation changes. A shared integer-sample offset across all six axes is estimated independently of EXG by averaging normalized cross-correlations over the central LSL excerpt. The selected offset is 55.102 seconds in nominal sample-index time; it is not an absolute timestamp correction and should not be applied to EXG or other recordings.

Each source and axis is standardized separately over the displayed window. LSL header units are g for acceleration and deg/s for angular velocity; decoded SD values are raw counts without established physical calibration. Consequently, these plots compare waveform shape, not physical amplitude accuracy. Residuals are standardized LSL minus standardized SD, with a common vertical scale across axes. Filtering, remaining timing differences, and noise can contribute to the residual, as discussed below for EXG.

For this excerpt, Pearson r ranges from 0.99930 to 0.99990 across the six axes. Alignment and evaluation use overlapping data, so this is a descriptive comparison, not independent validation. The bundled `data/imu_raw_sample.npz` retains original values for LSL 96–114 seconds and SD 145–180 seconds from the same source XDF files used for EXG. `data/imu_provenance.json` records hashes, stream names, sample ranges, labels, and source units. The example runs without the full recordings. `results/imu_comparison_arrays.npz` contains aligned standardized samples and residuals; `results/imu_summary.json` records the alignment and aggregate metrics.

## Spectral analysis: EXG and IMU

Run spectral analysis on the saved aligned, standardized 100–110 second segments:

```sh
# Refresh the aligned arrays if you change comparison processing:
.venv/bin/python compare.py
.venv/bin/python compare_imu.py
.venv/bin/python spectral_analysis.py
```

Use `--signal exg` or `--signal imu` to analyze only one signal type. The default is both.

| Analysis | EXG | IMU |
| --- | --- | --- |
| LSL and SD power spectral density | [PSD](results/exg_psd_10s.svg) | [PSD](results/imu_psd_10s.svg) |
| Residual power spectral density | [Residual PSD](results/exg_residual_psd_10s.svg) | [Residual PSD](results/imu_residual_psd_10s.svg) |
| Magnitude-squared coherence | [Coherence](results/exg_coherence_10s.svg) | [Coherence](results/imu_coherence_10s.svg) |

PNG versions are also generated. Power spectral density (PSD) shows how signal power is distributed across frequencies. Similar LSL and SD curves indicate similar spectral shapes. The separate residual PSD is computed from **LSL z-score minus SD z-score**, not by subtracting their power spectra; it shows which frequencies contribute to the remaining disagreement.

Coherence measures the consistency of the linear relationship between the recordings at each frequency, on a scale from 0 to 1. A value near 1 does not establish equal amplitudes, zero delay, or correct calibration. Interpret it alongside PSD: estimates where signal power is very small are less informative. Zero-power bins have undefined coherence and are saved as NaN.

The implementation uses [SciPy Welch spectral estimation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html), with two-second Hann windows, 50% overlap, mean removal within each window, and nine averaged segments. At 500 Hz this gives 0.5 Hz frequency-bin spacing; the window limits the ability to resolve nearby components. The short, overlapping segments provide a descriptive estimate, not a statistical significance test.

All PSD units are **z-score²/Hz**, because each displayed trace was standardized independently. These are comparisons of relative spectral content, not µV²/Hz, g²/Hz, or physical noise floors. EXG already has the 20 Hz high-pass filter applied; its shaded region below 20 Hz is attenuated by preprocessing. IMU has no high-pass filter and uses a logarithmic frequency axis to show slower motion. Both cover positive frequencies through the 250 Hz Nyquist limit; DC is retained in the saved arrays but omitted from plots. Anti-alias filtering affects the upper frequencies, and frequencies above 250 Hz in the original SD recordings are not analyzed here.

`results/exg_spectral_arrays.npz` and `results/imu_spectral_arrays.npz` store frequency bins, both PSDs, residual PSD, and coherence. The corresponding `*_spectral_metadata.json` files record the input array file, window settings, interval, units, and preprocessing. Spectral estimates reuse the existing alignment; they do not independently validate the decoder or identify the cause of a discrepancy.

## Included data

`data/raw_sample.npz` contains only 18 seconds of LSL data (96–114 s, 500 Hz, µV) and 20 seconds of SD data (150–170 s, 4,000 Hz, raw counts). Padding supports filtering and alignment of the central 10-second example. The arrays are `lsl` and `sd`, both samples × 16 channels. Source filenames, SHA-256 hashes, sample ranges, units, and rates are recorded in `data/provenance.json`. Times refer to sample indices divided by nominal rates, not synchronized absolute XDF timestamps.

## Processing and interpretation

1. Anti-alias FIR decimation reduces SD from 4,000 Hz to 500 Hz using zero-phase decimation by 8.
2. Both signals receive a fourth-order 20 Hz Butterworth high-pass applied forward and backward. Padding reduces boundary effects.
3. Normalized cross-correlation searches the excerpt for one shared offset across channels. It uses an interior LSL segment, excluding two seconds at each edge.
4. The aligned 100–110 s interval is extracted. Each channel and recording is standardized independently to zero mean and unit population standard deviation over those 10 seconds.
5. Pearson r, r², standardized RMSE, and a linear fit from LSL µV to SD counts are calculated per channel. The linear fit's RMSE is in counts. Fitted gains are empirical comparisons, not physical calibration.

Standardization supports waveform comparison but removes absolute amplitude differences; consult the fitted gains and residuals in the CSV for that aspect. Alignment and scores use overlapping data, so these are descriptive scores, not independent validation. This excerpt lies within one stable-offset segment. The full recording has several timing steps and needs local alignment; one offset must not be assumed for its entire duration. Cropped filtering may produce tiny numerical differences from filtering the complete files.

### What causes the residual error?

The [residual plot](results/residual_all_channels_10s.svg) shows **standardized LSL minus standardized SD** after filtering and alignment. Zero means the standardized traces agree at that sample. These differences are in z-score units, not µV or raw ADC counts; they measure waveform agreement after each trace's mean and amplitude scale have been removed.

**Downsampling can contribute, but this comparison does not isolate or quantify its contribution.** Reducing SD from 4,000 Hz to 500 Hz applies an anti-alias FIR filter before decimation, attenuating high-frequency content. If the LSL acquisition path uses a different filter response, the resulting waveforms can differ. Downsampling does not necessarily cause a mismatch when both paths have compatible bandwidth and timing. The LSL acquisition filter response has not been established here.

Other possible contributors include a remaining fractional-sample timing offset, clock drift within the window, differences in acquisition filtering, and noise. The current alignment uses one shared integer-sample shift; at 500 Hz, one sample is 2 ms. Even a smaller offset can produce visible residuals around sharp peaks. Independent standardization also affects the displayed differences, so this plot does not measure absolute amplitude accuracy.

To assess the downsampling contribution, first establish the acquisition filter responses, compare the signals over a matched bandwidth, and check how fractional-sample alignment and local clock correction change the residual. Those checks have not been performed in this example. **The residual alone does not establish a decoder error or show that downsampling is its main cause.**

## Files

- `decode_xf2.py`: primary XF2-to-XDF/EDF decoder, round-trip checks, and trace exports.
- `compare.py`: comparison processing, metrics, and plotting.
- `data/raw_sample.npz`, `data/provenance.json`: minimal raw example and source provenance.
- `results/aligned_all_channels_10s.png` / `.svg`: all channels with scores.
- `results/residual_all_channels_10s.png` / `.svg`: per-channel residual errors on a common z-score scale.
- `results/channel_metrics.csv`, `results/summary.json`: numerical results and alignment.
- `results/convolution_crosscorrelation.png`: curves for each channel, zoomed to ±100 ms.
- `results/comparison_arrays.npz`: aligned standardized traces, per-sample `residual_z` (LSL minus SD), and full convolution/cross-correlation arrays.

The curve labeled **Time-reversed overlap (convolution)** describes multiplication and summation with one signal reversed in time. **Time-shift similarity (cross-correlation)** measures agreement at different shifts. Convolution and cross-correlation operate on standardized 10-second signals and are divided by N. Nonzero-lag correlation uses fixed-N normalization, not overlap-specific Pearson normalization. The convolution plot centers its axis for display; the saved array includes the true output time axis. Convolution is not itself a similarity score.
