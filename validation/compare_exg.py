"""Reproduce the 10-second, 16-channel comparison using only bundled raw samples."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
from scipy import signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alpha', type=float, default=.45, help='SD opacity: 0 transparent, 1 opaque (default .45)')
    args = parser.parse_args()
    if not 0 < args.alpha <= 1:
        parser.error('--alpha must be in (0, 1]')
    out = ROOT / 'results'
    out.mkdir(exist_ok=True)
    meta = json.loads((ROOT / 'data/provenance.json').read_text())
    low, high = meta['sources']
    fs = low['sample_rate_hz']
    ratio = int(high['sample_rate_hz'] / fs)
    with np.load(ROOT / 'data/raw_sample.npz', allow_pickle=False) as data:
        x = data['lsl'].astype(float)
        y = signal.decimate(data['sd'].astype(float), ratio, ftype='fir', zero_phase=True, axis=0)
    assert x.shape[1] == y.shape[1] == 16
    assert np.isfinite(x).all() and np.isfinite(y).all()
    sos = signal.butter(4, 20, fs=fs, btype='highpass', output='sos')
    x = signal.sosfiltfilt(sos, x, axis=0)
    y = signal.sosfiltfilt(sos, y, axis=0)
    # Estimate shared lag on a central segment; preserve filter context at both ends.
    trim = int(2 * fs)
    template = x[trim:-trim]
    score = np.zeros(len(y) - len(template) + 1)
    for c in range(16):
        u = template[:, c] - template[:, c].mean()
        v = y[:, c]
        ones = np.ones(len(u))
        sm = signal.correlate(v, ones, mode='valid', method='fft')
        sq = signal.correlate(v*v, ones, mode='valid', method='fft') - sm*sm/len(u)
        score += signal.correlate(v, u, mode='valid', method='fft') / np.sqrt(np.maximum(sq, 1e-20) * np.sum(u*u))
    lag = int(score.argmax()) - trim
    start_s, end_s = meta['display_lsl_seconds']
    begin = int(start_s * fs) - low['start_sample']
    n = int((end_s - start_s) * fs)
    assert 0 <= begin+lag and begin+lag+n <= len(y)
    a, b = x[begin:begin+n], y[begin+lag:begin+lag+n]
    az = (a-a.mean(0))/a.std(0)
    bz = (b-b.mean(0))/b.std(0)
    t = start_s + np.arange(n)/fs
    rows = []
    for c in range(16):
        r = np.corrcoef(a[:, c], b[:, c])[0, 1]
        gain = np.dot(a[:, c]-a[:, c].mean(), b[:, c]-b[:, c].mean()) / np.sum((a[:, c]-a[:, c].mean())**2)
        intercept = b[:, c].mean() - gain*a[:, c].mean()
        rows.append(dict(lsl_channel=f'Ch{c}', sd_channel=f'Ch {c+1:02d}', pearson_r=float(r), r_squared=float(r*r), z_score_rmse=float(np.sqrt(np.mean((az[:, c]-bz[:, c])**2))), gain_counts_per_uV=float(gain), fit_rmse_counts=float(np.sqrt(np.mean((b[:, c]-gain*a[:, c]-intercept)**2)))))
    with (out/'channel_metrics.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    offset = high['start_sample']/ratio - low['start_sample'] + lag
    summary = dict(display_lsl_seconds=[start_s, end_s],samples_per_channel=n, sample_rate_hz=fs, global_sd_minus_lsl_samples=offset, global_sd_minus_lsl_seconds=offset/fs, pearson_r_min=min(r['pearson_r'] for r in rows), pearson_r_max=max(r['pearson_r'] for r in rows), sd_alpha=args.alpha)
    (out/'summary.json').write_text(json.dumps(summary, indent=2))
    residual = az - bz
    residual_limit = max(0.1, float(np.max(np.abs(residual))) * 1.05)
    fig, axes = plt.subplots(16, 1, figsize=(20, 38), sharex=True, layout='constrained')
    for c, ax in enumerate(axes):
        ax.plot(t, az[:, c], color='#1769aa', lw=1.1, alpha=.95, label='LSL')
        ax.plot(t, bz[:, c], color='#f07818', lw=.9, alpha=args.alpha, label='SD, downsampled and aligned')
        ax.set_title(f"LSL Ch{c} ↔ SD Ch {c+1:02d}     |     Pearson r = {rows[c]['pearson_r']:.4f}     |     RMSE (z-score) = {rows[c]['z_score_rmse']:.3f}", loc='left', fontsize=14, pad=6)
        ax.set_ylabel('Amplitude\n(z-score)', fontsize=12, labelpad=12)
        ax.tick_params(labelsize=11)
        ax.grid(alpha=.18)
        ax.set_xlim(start_s, end_s)
    axes[0].legend(loc='upper right', fontsize=12)
    axes[-1].set_xlabel('LSL sample-index time (s)', fontsize=15)
    fig.suptitle('All 16 channels · aligned 100–110 s · 20 Hz high-pass\nPer-channel scores computed on these 10 seconds; each trace standardized separately', fontsize=20)
    for ext in ['png','svg']:
        fig.savefig(out/f'aligned_all_channels_10s.{ext}', dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(16, 1, figsize=(20, 38), sharex=True, sharey=True,
                             layout='constrained')
    for c, ax in enumerate(axes):
        ax.axhline(0, color='#555555', lw=.8)
        ax.plot(t, residual[:, c], color='#7b3294', lw=.8)
        ax.set_title(f"LSL Ch{c} − SD Ch {c+1:02d}     |     RMSE (z-score) = {rows[c]['z_score_rmse']:.3f}",
                     loc='left', fontsize=14, pad=6)
        ax.set_ylabel('Residual\n(z-score)', fontsize=12, labelpad=12)
        ax.set_ylim(-residual_limit, residual_limit)
        ax.set_xlim(start_s, end_s)
        ax.tick_params(labelsize=11)
        ax.grid(alpha=.18)
    axes[-1].set_xlabel('LSL sample-index time (s)', fontsize=15)
    fig.suptitle('Residual error · aligned 100–110 s · 20 Hz high-pass\n'
                 'Standardized LSL − standardized SD; zero indicates agreement\n'
                 'Common vertical scale across channels; units are z-scores, not µV or ADC counts',
                 fontsize=20)
    for ext in ['png', 'svg']:
        fig.savefig(out/f'residual_all_channels_10s.{ext}', dpi=170)
    plt.close(fig)
    lags = signal.correlation_lags(n,n)/fs
    corr = np.column_stack([signal.correlate(bz[:, c],az[:, c],method='fft')/n for c in range(16)])
    conv = np.column_stack([signal.fftconvolve(az[:, c],bz[:, c])/n for c in range(16)])
    np.savez_compressed(out/'comparison_arrays.npz',time_s=t,lsl_z=az,sd_z=bz,crosscorrelation=corr,lag_s=lags,residual_z=residual,convolution=conv,convolution_time_s=np.arange(2*n-1)/fs)
    fig, axes = plt.subplots(4,4,figsize=(16,12),layout='constrained')
    for c,ax in enumerate(axes.flat):
        ax.plot(lags,corr[:,c],label='Time-shift similarity (cross-correlation)')
        ax.plot(lags,conv[:,c],alpha=.5,label='Time-reversed overlap (convolution)')
        ax.set(title=f'Ch{c}',xlabel='Shift / centered reversed-overlap time (s)',xlim=(-.1,.1))
    axes.flat[0].legend(fontsize=8)
    fig.suptitle('10-second sample · standardized signals · sums divided by N\nTime-reversed overlap compares one signal with the other reversed; it is not an alignment score')
    fig.savefig(out/'convolution_crosscorrelation.png',dpi=150)
    plt.close(fig)
    print(json.dumps(summary,indent=2))

if __name__ == '__main__':
    main()
