"""End-to-end example: decode an XF2 file with decode_xf2, then read the exports.

Usage:
    python examples/quickstart.py                    # builds a small synthetic XF2
    python examples/quickstart.py recording.xf2      # uses your own recording

Outputs go to examples/output/ (ignored by git).
"""
from pathlib import Path
import json
import struct
import subprocess
import sys
import numpy as np
import pyedflib
import pyxdf

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import decode_xf2  # noqa: E402


def synthetic_xf2(path, seconds=2):
    """Write EXG (2 ch, 500 Hz) and IMU (6 axes, 100 Hz) records in the observed layout."""
    blob = bytearray()
    for index in range(10 * seconds):
        stamp = 1700000000 + index / 10
        sec, ms = int(stamp), round((stamp % 1) * 1000)
        t = stamp + np.arange(50) / 500
        exg = 32768 + 1000 * np.column_stack([np.sin(2*np.pi*10*t), np.cos(2*np.pi*25*t)])
        imu = np.tile([0, 0, 16384, 10, -10, 0], (10, 1))
        for kind, mapping, rate, values, dtype in ((160, 0b11, 500, exg, '<u2'),
                                                   (161, 0b11, 100, imu, '>i2')):
            payload = np.asarray(values).astype(dtype).tobytes()
            blob += b'\r' + struct.pack('<BIHH', kind, sec, ms, 6 + len(payload))
            blob += struct.pack('<HHHB', index, mapping, rate, 1)
            blob += payload + b'\0\0\n'  # two CRC bytes (unverified), then terminator
    path.write_bytes(blob)


def main():
    out = ROOT / 'examples' / 'output'
    out.mkdir(parents=True, exist_ok=True)
    if len(sys.argv) > 1:
        source = Path(sys.argv[1]).expanduser()
    else:
        source = out / 'synthetic.xf2'
        synthetic_xf2(source)

    # 1. Python API: decode into in-memory arrays of raw counts.
    streams, report = decode_xf2.decode(source)
    for s in streams:
        print(f"{s['name']:>18}: {s['values'].shape[0]} samples x {len(s['labels'])} "
              f"channels at {s['rate']} Hz nominal")
    print('Clock fit:', json.dumps(report['streams'], indent=2))

    # 2. Command line: export verified XDF and EDF+ files.
    subprocess.run([sys.executable, str(ROOT / 'decode_xf2.py'), str(source),
                    '--format', 'both', '--output-dir', str(out / 'decoded')],
                   check=True, stdout=subprocess.DEVNULL)

    # 3. Read the exports back with standard tools.
    xdf, _ = pyxdf.load_xdf(str(out / 'decoded' / (source.stem + '.xdf')),
                            synchronize_clocks=False, dejitter_timestamps=False)
    for s in xdf:
        print('XDF stream', s['info']['name'][0], s['time_series'].shape)
    exg_edf = out / 'decoded' / (source.stem + '_EXG.edf')
    if exg_edf.exists():
        with pyedflib.EdfReader(str(exg_edf)) as reader:
            print('EDF EXG labels', reader.getSignalLabels(),
                  'first counts', reader.readSignal(0)[:5])
    print('Outputs written to', out / 'decoded')


if __name__ == '__main__':
    main()
