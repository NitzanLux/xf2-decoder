"""Compare bundled LSL and SD IMU excerpts without removing low-frequency motion."""
from pathlib import Path
import csv
import json
import numpy as np
from scipy import signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent


def main():
    out = ROOT / 'results' / 'imu'
    out.mkdir(parents=True, exist_ok=True)
    meta = json.loads((ROOT / 'data/imu_provenance.json').read_text())
    low, high = meta['sources']
    fs = low['sample_rate_hz']
    ratio = high['sample_rate_hz'] / fs
    if ratio != int(ratio) or ratio < 1:
        raise ValueError('SD rate must be an integer multiple of LSL rate')
    with np.load(ROOT / 'data/imu_raw_sample.npz', allow_pickle=False) as data:
        x, raw_sd = data['lsl'].astype(float), data['sd'].astype(float)
    if not np.isfinite(x).all() or not np.isfinite(raw_sd).all():
        raise ValueError('IMU excerpts contain missing/nonfinite samples; explicit gap handling required')
    y = signal.decimate(raw_sd, int(ratio), ftype='fir', zero_phase=True, axis=0)
    labels = [c['label'] for c in low['channels']]
    assert labels == [c['label'] for c in high['channels']]
    assert x.shape[1] == y.shape[1] == len(labels)
    # One shared integer shift across axes; center each candidate window.
    trim = int(2 * fs)
    template = x[trim:-trim]
    scores = np.zeros(len(y) - len(template) + 1)
    active = 0
    for c in range(len(labels)):
        u = template[:, c] - template[:, c].mean()
        energy = np.dot(u, u)
        if energy <= 1e-20:
            continue
        v = y[:, c]
        ones = np.ones(len(u))
        sums = signal.correlate(v, ones, mode='valid', method='fft')
        variance = signal.correlate(v*v, ones, mode='valid', method='fft') - sums*sums/len(u)
        scores += signal.correlate(v, u, mode='valid', method='fft') / np.sqrt(np.maximum(variance, 1e-20)*energy)
        active += 1
    if not active:
        raise ValueError('No varying IMU axes available for alignment')
    peak = int(scores.argmax())
    if peak in (0, len(scores)-1):
        raise ValueError('Alignment is at search boundary; supply a wider SD excerpt')
    lag = peak - trim
    start, stop = meta['display_lsl_seconds']
    begin = int(start*fs) - low['start_sample']
    n = int((stop-start)*fs)
    assert 0 <= begin+lag and begin+lag+n <= len(y)
    a, b = x[begin:begin+n], y[begin+lag:begin+lag+n]
    if np.any(a.std(0) == 0) or np.any(b.std(0) == 0):
        raise ValueError('Cannot standardize a constant IMU axis')
    az, bz = (a-a.mean(0))/a.std(0), (b-b.mean(0))/b.std(0)
    residual = az-bz
    t = start + np.arange(n)/fs
    rows = [dict(channel=label, pearson_r=float(np.mean(az[:, c]*bz[:, c])),
                 z_score_rmse=float(np.sqrt(np.mean(residual[:, c]**2))))
            for c, label in enumerate(labels)]
    with (out/'channel_metrics.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for errors in (False, True):
        fig, axes = plt.subplots(len(labels), 1, figsize=(18, 16), sharex=True,
                                 sharey=errors, layout='constrained')
        for c, ax in enumerate(axes):
            if errors:
                ax.axhline(0, color='#555555', lw=.8)
                ax.plot(t, residual[:, c], color='#7b3294', lw=.8)
            else:
                ax.plot(t, az[:, c], label='LSL', color='#1769aa', lw=1.1)
                ax.plot(t, bz[:, c], label='SD, downsampled and aligned', color='#f07818', alpha=.45, lw=.9)
            ax.set_title(f"{labels[c]} | Pearson r = {rows[c]['pearson_r']:.4f} | RMSE (z-score) = {rows[c]['z_score_rmse']:.3f}", loc='left')
            ax.set_ylabel('Residual (z-score)' if errors else 'Amplitude (z-score)')
            ax.set_xlim(start, stop)
            ax.grid(alpha=.18)
        if not errors:
            axes[0].legend(loc='upper right')
        axes[-1].set_xlabel('LSL IMU sample-index time (s)')
        title = 'Residual: standardized LSL − standardized SD; common vertical scale' if errors else 'Aligned traces; each axis and source standardized independently'
        fig.suptitle(f'IMU · {start}–{stop} s · {title}\nSD anti-alias decimation only; no 20 Hz high-pass; not a physical calibration comparison')
        stem = 'residual_all_channels_10s' if errors else 'aligned_all_channels_10s'
        for ext in ('png', 'svg'):
            fig.savefig(out/f'{stem}.{ext}', dpi=150)
        plt.close(fig)
    np.savez_compressed(out/'comparison_arrays.npz', time_s=t, lsl_z=az, sd_z=bz,
                        residual_z=residual, alignment_scores=scores/active)
    summary = dict(display_lsl_seconds=[start, stop], samples_per_channel=n,
                   sample_rate_hz=fs, global_sd_minus_lsl_seconds=(high['start_sample']/ratio-low['start_sample']+lag)/fs,
                   alignment_peak_mean_correlation=float(scores[peak]/active),
                   pearson_r_min=min(r['pearson_r'] for r in rows),
                   pearson_r_max=max(r['pearson_r'] for r in rows),
                   note='Independent IMU integer-sample alignment; descriptive same-window comparison, not calibration or absolute clock synchronization.')
    (out/'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
