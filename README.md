# XDF channel comparison — standalone example

A reproducible comparison of 16 LSL channels and their corresponding SD channels over **LSL sample-index time 100–110 seconds**. Includes small original-value recording excerpts, analysis code, per-channel metrics, and generated figures. No full recordings or external data paths are required.

![All aligned channels](results/aligned_all_channels_10s.png)

Open [the SVG](results/aligned_all_channels_10s.svg) for lossless zoom, or [the CSV](results/channel_metrics.csv) for quantitative results. Each channel title reports Pearson r and standardized RMSE for the actual displayed 10 seconds.

## Run

Python 3.9 or later:

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
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

- `compare.py`: all processing, metrics, and plotting.
- `data/raw_sample.npz`, `data/provenance.json`: minimal raw example and source provenance.
- `results/aligned_all_channels_10s.png` / `.svg`: all channels with scores.
- `results/channel_metrics.csv`, `results/summary.json`: numerical results and alignment.
- `results/convolution_crosscorrelation.png`: curves for each channel, zoomed to ±100 ms.
- `results/comparison_arrays.npz`: aligned standardized traces and full convolution/cross-correlation arrays.

Convolution and cross-correlation operate on standardized 10-second signals and are divided by N. Nonzero-lag correlation uses fixed-N normalization, not overlap-specific Pearson normalization. The convolution plot centers its axis for display; the saved array includes the true output time axis. Convolution is not itself a similarity score.
