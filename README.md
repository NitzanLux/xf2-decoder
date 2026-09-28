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

Open [the SVG](results/aligned_all_channels_10s.svg) for lossless zoom, or [the CSV](results/channel_metrics.csv) for quantitative results. Each channel title reports Pearson r and standardized RMSE for the actual displayed 10 seconds.

## Run the comparison

After installing the dependencies:

```sh
.venv/bin/python compare.py
```

The SD curve uses alpha 0.45 and LSL uses 0.95, revealing both overlapping traces. Alpha is **opacity**: lower values give more transparency. To increase SD opacity:

```sh
.venv/bin/python compare.py --alpha 0.9
```

## Included data

`data/raw_sample.npz` contains only 18 seconds of LSL data (96–114 s, 500 Hz, µV) and 20 seconds of SD data (150–170 s, 4,000 Hz, raw counts). Padding supports filtering and alignment of the central 10-second example. The arrays are `lsl` and `sd`, both samples × 16 channels. Source filenames, SHA-256 hashes, sample ranges, units, and rates are recorded in `data/provenance.json`. Times refer to sample indices divided by nominal rates, not synchronized absolute XDF timestamps.

## Processing and interpretation

1. Anti-alias FIR decimation reduces SD from 4,000 Hz to 500 Hz using zero-phase decimation by 8.
2. Both signals receive a fourth-order 20 Hz Butterworth high-pass applied forward and backward. Padding reduces boundary effects.
3. Normalized cross-correlation searches the excerpt for one shared offset across channels. It uses an interior LSL segment, excluding two seconds at each edge.
4. The aligned 100–110 s interval is extracted. Each channel and recording is standardized independently to zero mean and unit population standard deviation over those 10 seconds.
5. Pearson r, r², standardized RMSE, and a linear fit from LSL µV to SD counts are calculated per channel. The linear fit's RMSE is in counts. Fitted gains are empirical comparisons, not physical calibration.

Standardization supports waveform comparison but removes absolute amplitude differences; consult the fitted gains and residuals in the CSV for that aspect. Alignment and scores use overlapping data, so these are descriptive scores, not independent validation. This excerpt lies within one stable-offset segment. The full recording has several timing steps and needs local alignment; one offset must not be assumed for its entire duration. Cropped filtering may produce tiny numerical differences from filtering the complete files.

## Files

- `decode_xf2.py`: primary XF2-to-XDF/EDF decoder, round-trip checks, and trace exports.
- `compare.py`: comparison processing, metrics, and plotting.
- `data/raw_sample.npz`, `data/provenance.json`: minimal raw example and source provenance.
- `results/aligned_all_channels_10s.png` / `.svg`: all channels with scores.
- `results/channel_metrics.csv`, `results/summary.json`: numerical results and alignment.
- `results/convolution_crosscorrelation.png`: curves for each channel, zoomed to ±100 ms.
- `results/comparison_arrays.npz`: aligned standardized traces and full convolution/cross-correlation arrays.

Convolution and cross-correlation operate on standardized 10-second signals and are divided by N. Nonzero-lag correlation uses fixed-N normalization, not overlap-specific Pearson normalization. The convolution plot centers its axis for display; the saved array includes the true output time axis. Convolution is not itself a similarity score.
