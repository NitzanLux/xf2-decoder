"""Decode the observed XtrodES XF2 layout; preserve ADC counts and packet timing.

Usage: python decode_xf2.py "distance check"
Dependencies: numpy matplotlib pyxdf
"""
from pathlib import Path
import argparse
import collections
import hashlib
import json
import struct
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyxdf


def decode(path):
    blob = path.read_bytes()
    pos = padding = 0
    records = collections.defaultdict(list)
    events = []
    while pos < len(blob):
        if blob[pos] == 0:
            pos += 1
            padding += 1
            continue
        assert blob[pos] == 13, f'Bad start at {pos}'
        kind, sec, ms, length = struct.unpack_from('<BIHH', blob, pos + 1)
        assert ms < 1000
        assert kind in (0, 2, 160, 161), f'Unsupported record {kind} at {pos}'
        sensor = kind in (160, 161)
        end = pos + 13 + length + int(sensor)
        assert end <= len(blob) and blob[end - 1] == 10, f'Bad end at {pos}'
        if sensor:
            index, mapping, rate, down = struct.unpack_from('<HHHB', blob, pos + 10)
            channels = bin(mapping).count("1") * (3 if kind == 161 else 1)
            assert down > 0 and (length - 6) % (2 * channels) == 0
            # XF2 length counts the six configuration bytes, excluding downsample.
            raw = np.frombuffer(blob[pos + 17:end - 3], dtype='>i2' if kind == 161 else '<u2')
            values = raw.reshape(-1, channels).astype(np.int32)
            records[kind].append((sec + ms / 1000, index, mapping, rate / down, values))
        else:
            events.append(dict(type=kind, device_unix_seconds=sec + ms / 1000,
                               offset=pos, payload_hex=blob[pos + 10:end - 3].hex()))
        pos = end
    streams = []
    stats = {}
    for kind, recs in sorted(records.items()):
        name = 'EXG' if kind == 160 else 'IMU'
        assert len({(r[2], r[3], r[4].shape) for r in recs}) == 1
        anchors = np.array([r[0] for r in recs])
        indices = np.array([r[1] for r in recs])
        gaps = (np.diff(indices) % 65536) - 1
        assert np.all(gaps == 0), 'Packet loss requires explicit gap handling'
        values = np.concatenate([r[4] for r in recs])
        count = len(recs[0][4])
        x = np.arange(len(recs)) * count
        slope, intercept = np.polyfit(x, anchors - anchors[0], 1)
        stamps = anchors[0] + intercept + np.arange(len(values)) * slope
        residual = anchors - (anchors[0] + intercept + x * slope)
        labels = [f'Ch {i+1:02d}' for i in range(values.shape[1])] if kind == 160 else ['Acc X', 'Acc Y', 'Acc Z', 'Gyro X', 'Gyro Y', 'Gyro Z']
        assert kind != 161 or recs[0][2] == 3
        stats[name] = dict(samples=len(values), channels=values.shape[1], packets=len(recs),
                           nominal_rate=recs[0][3], fitted_rate=1/slope,
                           duration_seconds=float(stamps[-1]-stamps[0]),
                           first_packet_index=int(indices[0]), last_packet_index=int(indices[-1]),
                           missing_packets=int(gaps.sum()), first_device_timestamp=float(anchors[0]),
                           last_device_timestamp=float(anchors[-1]),
                           max_clock_fit_residual_ms=float(np.max(np.abs(residual))*1000))
        streams.append(dict(name=name, values=values, stamps=stamps, labels=labels,
                            rate=recs[0][3], fmt='int32', unit='count'))
        packet_values = np.column_stack((indices, x, np.full(len(recs), count)))
        streams.append(dict(name=name+'_packet_timing', values=packet_values.astype(np.int32),
                            stamps=anchors, labels=['packet_index', 'first_sample_index', 'sample_count'],
                            rate=0, fmt='int32', unit='count'))
    return streams, dict(source=path.name, sha256=hashlib.sha256(blob).hexdigest(),
                        padding_bytes=padding, events=events, streams=stats,
                        crc='Stored CRC bytes present; algorithm not verified')


def varint(value):
    return b'\x04' + struct.pack('<I', value)


def chunk(f, tag, body):
    f.write(varint(len(body)+2) + struct.pack('<H', tag) + body)


def write_xdf(path, streams, report):
    with path.open('wb') as f:
        f.write(b'XDF:')
        chunk(f, 1, b'<info><version>1.0</version></info>')
        for sid, s in enumerate(streams, 1):
            root = ET.Element('info')
            for key, val in dict(name=s['name'], type=s['name'], channel_count=len(s['labels']),
                                 nominal_srate=s['rate'], channel_format=s['fmt'],
                                 source_id=report['sha256']+':'+s['name']).items():
                ET.SubElement(root, key).text = str(val)
            desc = ET.SubElement(root, 'desc')
            ET.SubElement(desc, 'source_file').text = report['source']
            ET.SubElement(desc, 'timestamp_basis').text = 'Device Unix seconds; not synchronized to LSL host clock. Sample timestamps fit packet anchors by sample index; original anchors retained in packet_timing streams. Anchor treated as first sample; device latency unknown.'
            ET.SubElement(desc, 'calibration').text = 'Unscaled raw counts. EXG unsigned 16-bit little-endian; IMU signed 16-bit big-endian. Physical calibration not established.'
            ET.SubElement(desc, 'conversion_report').text = json.dumps(report)
            channels = ET.SubElement(desc, 'channels')
            for label in s['labels']:
                channel = ET.SubElement(channels, 'channel')
                ET.SubElement(channel, 'label').text = label
                ET.SubElement(channel, 'unit').text = s['unit']
            chunk(f, 2, struct.pack('<I', sid)+ET.tostring(root))
        for sid, s in enumerate(streams, 1):
            dtype = np.dtype([('flag', 'u1'), ('timestamp', '<f8'), ('values', '<i4', (len(s['labels']),))])
            for start in range(0, len(s['values']), 20000):
                values = s['values'][start:start+20000]
                packed = np.empty(len(values), dtype=dtype)
                packed['flag'] = 8
                packed['timestamp'] = s['stamps'][start:start+len(values)]
                packed['values'] = values
                chunk(f, 3, struct.pack('<I', sid)+varint(len(values))+packed.tobytes())
            footer = f'<info><first_timestamp>{s["stamps"][0]}</first_timestamp><last_timestamp>{s["stamps"][-1]}</last_timestamp><sample_count>{len(s["values"])}</sample_count></info>'
            chunk(f, 6, struct.pack('<I', sid)+footer.encode())


def envelope(t, y, bins=2200):
    step = max(1, len(y)//bins)
    starts = np.arange(0, len(y), step)
    ends = np.minimum(starts + step, len(y)) - 1
    tt = (t[starts] + t[ends]) / 2
    return tt, np.minimum.reduceat(y, starts), np.maximum.reduceat(y, starts)


def plot_stream(path, stem, s):
    n = len(s['labels'])
    fig, axes = plt.subplots(n, 1, figsize=(15, 1.25*n+1.6), sharex=True)
    t = s['stamps'] - s['stamps'][0]
    for i, ax in enumerate(axes):
        tt, low, high = envelope(t, s['values'][:, i])
        ax.fill_between(tt, low, high, color='#156b9a', alpha=.85, linewidth=.3)
        ax.set_ylabel(s['labels'][i], rotation=0, ha='right', fontsize=9)
        ax.grid(alpha=.18)
        ax.tick_params(labelsize=8)
    axes[-1].set_xlabel('Elapsed seconds from stream start')
    fig.suptitle(f'{stem} | {s["name"]} traces | raw counts\nUnfiltered min/max envelope; independent channel scales', fontsize=14)
    fig.tight_layout(rect=(0,0,1,.965))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('folder', type=Path)
    args = parser.parse_args()
    out = args.folder / 'decoded'
    out.mkdir(exist_ok=True)
    reports = []
    for source in sorted(args.folder.glob('*.xf2')):
        streams, report = decode(source)
        target = out / (source.stem+'.xdf')
        write_xdf(target, streams, report)
        loaded, _ = pyxdf.load_xdf(str(target), synchronize_clocks=False, dejitter_timestamps=False)
        assert len(loaded) == len(streams)
        for actual, expected in zip(loaded, streams):
            assert np.array_equal(actual['time_series'], expected['values'])
            assert np.array_equal(actual['time_stamps'], expected['stamps'])
            assert np.all(np.diff(actual['time_stamps']) > 0)
        report['xdf_roundtrip_verified'] = True
        for s in streams:
            if s['rate']:
                plot_stream(out/(source.stem+'_'+s['name']+'_traces.png'), source.stem, s)
        reports.append(report)
        print(json.dumps(report), flush=True)
    (out/'conversion_report.json').write_text(json.dumps(reports, indent=2))
    (out/'README.txt').write_text(
        'Each XDF contains EXG, IMU, and two packet-timing streams. Counts are raw, unfiltered, and uncalibrated.\n'
        'EXG stores unsigned 16-bit source values losslessly in int32. IMU uses signed big-endian int16 source values.\n'
        'EXG nominal rate: 4000 Hz; IMU: 1000 Hz. Sample clocks fitted independently to device packet timestamps.\n'
        'Timestamps are device Unix seconds, not synchronized LSL host time. Packet anchors assumed to indicate first sample.\n'
        'Original packet timestamps and indices remain in the packet_timing streams. CRC algorithm was not verified.\n'
        'All frame boundaries, payload sizes, packet continuity and complete consumption except zero padding checked.\n'
        'Every XDF was reloaded with pyxdf and every sample and timestamp compared exactly.\n'
        '0002 continues 0001 by packet index; 0000 is a separate recording.\n'
        'PNG plots show min/max envelopes in raw counts with independent channel scales; no peaks discarded by simple subsampling.\n'
        'For pyxdf, use synchronize_clocks=False, dejitter_timestamps=False to retain exported timestamps exactly.\n')


if __name__ == '__main__':
    main()
