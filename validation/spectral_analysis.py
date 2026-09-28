"""Welch spectra and coherence for the saved aligned EXG and IMU comparisons."""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy import signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent


def spectra(lsl, sd, fs):
    """Use multiple two-second Hann windows, 50% overlap, and mean removal."""
    lsl, sd = np.asarray(lsl), np.asarray(sd)
    if lsl.shape != sd.shape or lsl.ndim != 2:
        raise ValueError('Expected matching samples-by-channels arrays')
    if not np.isfinite(lsl).all() or not np.isfinite(sd).all() or fs <= 0:
        raise ValueError('Spectra require finite samples and a positive sample rate')
    nperseg = int(round(2 * fs))
    if nperseg < 2 or len(lsl) < 2*nperseg:
        raise ValueError('At least four seconds of samples required for averaged coherence')
    settings = dict(fs=fs, window='hann', nperseg=nperseg,
                    noverlap=nperseg//2, detrend='constant', axis=0)
    f, p_lsl = signal.welch(lsl, scaling='density', **settings)
    _, p_sd = signal.welch(sd, scaling='density', **settings)
    _, p_error = signal.welch(lsl-sd, scaling='density', **settings)
    _, cross = signal.csd(lsl, sd, scaling='density', **settings)
    denominator = p_lsl*p_sd
    coherence = np.full_like(denominator, np.nan)
    np.divide(np.abs(cross)**2, denominator, out=coherence, where=denominator > 0)
    coherence = np.clip(coherence, 0, 1)
    return dict(frequency_hz=f, lsl_psd=p_lsl, sd_psd=p_sd,
                residual_psd=p_error, coherence=coherence), dict(
                    sample_rate_hz=fs, samples=len(lsl), window='periodic Hann',
                    window_samples=nperseg, overlap_samples=nperseg//2,
                    averaged_segments=1+(len(lsl)-nperseg)//(nperseg-nperseg//2),
                    frequency_bin_spacing_hz=fs/nperseg,
                    detrend='Mean removed per window', psd_units='z-score squared / Hz')


def analyze(kind, out):
    prefix = 'imu_' if kind == 'imu' else ''
    source = out / f'{prefix}comparison_arrays.npz'
    with np.load(source, allow_pickle=False) as arrays:
        t, lsl, sd = arrays['time_s'], arrays['lsl_z'], arrays['sd_z']
    dt = np.diff(t)
    if not len(dt) or np.any(dt <= 0) or not np.allclose(dt, dt[0], atol=1e-10, rtol=1e-7):
        raise ValueError('Spectral input requires a uniform increasing time grid')
    fs = float(round(1 / np.median(dt), 8))
    values, metadata = spectra(lsl, sd, fs)
    labels = ['Acc X', 'Acc Y', 'Acc Z', 'Gyro X', 'Gyro Y', 'Gyro Z'] if kind == 'imu' else [f'LSL Ch{c} / SD Ch {c+1:02d}' for c in range(lsl.shape[1])]
    if len(labels) != lsl.shape[1]:
        raise ValueError('Unexpected number of IMU channels')
    metadata.update(source_arrays=source.name, channels=labels,
                    interval_seconds=[float(t[0]), float(t[-1]+1/fs)],
                    preprocessing='SD anti-alias decimation; separate per-channel standardization. ' +
                    ('No high-pass filter.' if kind == 'imu' else 'Both streams have a 20 Hz high-pass filter.'),
                    interpretation='Relative spectral shape, not physical calibration. Coherence is undefined at zero-power bins and uncertain where power is small.')
    np.savez_compressed(out/f'{kind}_spectral_arrays.npz', **values)
    (out/f'{kind}_spectral_metadata.json').write_text(json.dumps(metadata, indent=2))
    f = values['frequency_hz']
    positive = f > 0
    cols = 2 if kind == 'imu' else 4
    rows = (len(labels)+cols-1)//cols
    for measure in ('psd', 'residual_psd', 'coherence'):
        fig, axes = plt.subplots(rows, cols, figsize=(cols*4.6, rows*3),
                                 sharex=True, sharey=True, squeeze=False, layout='constrained')
        for c, ax in enumerate(axes.flat):
            if c >= len(labels):
                ax.set_visible(False)
                continue
            if measure == 'psd':
                ax.semilogy(f[positive], np.maximum(values['lsl_psd'][positive, c], 1e-20), label='LSL', color='#1769aa')
                ax.semilogy(f[positive], np.maximum(values['sd_psd'][positive, c], 1e-20), label='SD', color='#f07818', alpha=.75)
            elif measure == 'residual_psd':
                ax.semilogy(f[positive], np.maximum(values[measure][positive, c], 1e-20), color='#7b3294')
            else:
                ax.plot(f[positive], values[measure][positive, c], color='#247b52', lw=1)
                ax.set_ylim(0, 1.02)
            if kind == 'imu':
                ax.set_xscale('log')
            else:
                ax.axvspan(0, 20, color='grey', alpha=.12)
            ax.set_xlim(f[1], fs/2)
            ax.set_title(labels[c])
            ax.set_xlabel('Frequency (Hz)')
            ax.set_ylabel('Magnitude-squared coherence' if measure == 'coherence' else 'PSD (z-score²/Hz)')
            ax.grid(alpha=.2, which='both')
        if measure == 'psd':
            axes.flat[0].legend()
        title = dict(psd='LSL and SD power spectra', residual_psd='Residual power spectrum: standardized LSL − SD', coherence='LSL–SD magnitude-squared coherence')[measure]
        note = 'Log frequency axis; no high-pass filter' if kind == 'imu' else 'Shaded <20 Hz region attenuated by preprocessing'
        fig.suptitle(f'{kind.upper()} · {title} · aligned {t[0]:g}–{t[-1]+1/fs:g} s\n'
                     f'Welch: 2 s Hann windows, 50% overlap, {metadata["averaged_segments"]} segments\n{note}', fontsize=11)
        for ext in ('png', 'svg'):
            fig.savefig(out/f'{kind}_{measure}_10s.{ext}', dpi=150)
        plt.close(fig)
    print(f'{kind.upper()}: saved power spectra, separate residual spectra, and coherence')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--signal', choices=('exg', 'imu', 'both'), default='both')
    args = parser.parse_args()
    for kind in (('exg', 'imu') if args.signal == 'both' else (args.signal,)):
        analyze(kind, ROOT/'results')


if __name__ == '__main__':
    main()
