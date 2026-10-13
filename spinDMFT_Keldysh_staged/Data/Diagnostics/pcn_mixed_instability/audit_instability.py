#!/usr/bin/env python3
"""Summarize the retained one-rank mixed-stage instability controls."""
from pathlib import Path
import json
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent

def plain(value):
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value

records = []
for case in sorted(ROOT.iterdir()):
    if not case.is_dir():
        continue
    files = list(case.rglob('*.hdf5'))
    if len(files) != 1:
        continue
    with h5py.File(files[0]) as f:
        attrs = f['runtimedata'].attrs
        re = f['results/Re_correlation'][...]
        im = f['results/Im_correlation'][...]
        errors = np.hypot(f['runtimedata/Re_correlation_sample_stds'][...],
                          f['runtimedata/Im_correlation_sample_stds'][...])
        z = re + 1j * im
        row = {'case': case.name, 'file': str(files[0]),
               'real_time_started': plain(f['preparation'].attrs['real_time_started']),
               'preparation_termination': plain(f['preparation'].attrs['termination']),
               'termination': plain(attrs['termination']),
               'sampling_strategy': plain(f['parameters'].attrs['sampling_strategy']),
               'iterations': plain(attrs['num_Iterations']),
               'max_absolute_correlation': float(abs(z).max()),
               'max_correlation_error': float(errors.max()),
               'absolute_residuals': plain(attrs['absolute_iteration_errors']),
               'phase_magnitudes': plain(attrs['average_phase_magnitudes']),
               'closure_residuals': plain(attrs['denominator_constancy_residuals']),
               'times': []}
        times = f['results'].attrs['real_times']
        for target in [0, 1, 3, 5, 8, 10, 12, 15]:
            if target > times[-1]:
                continue
            i = int(np.argmin(abs(times - target)))
            row['times'].append({'time': float(times[i]),
                                 'max_absolute_correlation': float(abs(z[i]).max()),
                                 'max_correlation_error': float(errors[i].max())})
        records.append(row)
(ROOT / 'audit_results.json').write_text(json.dumps(records, indent=2) + '\n')
for row in records:
    if row['real_time_started'] != 'yes':
        print(row['case'], 'preparation not converged')
    else:
        print(row['case'], row['termination'], row['iterations'],
              'max|g|=', row['max_absolute_correlation'],
              'error=', row['max_correlation_error'],
              'phase=', row['phase_magnitudes'][-1])
