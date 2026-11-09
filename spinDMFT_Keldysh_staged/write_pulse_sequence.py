#!/usr/bin/env python3
"""Write HDF5 pulse times, axes and angles in radians (positive +i angle convention)."""
import argparse
from pathlib import Path

import h5py
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--pulse', nargs=3, action='append', default=[],
                        metavar=('TIME', 'AXIS', 'ANGLE'),
                        help='repeat in time order; AXIS is x, y, z, or e.g. 1,1,0')
    args = parser.parse_args()
    times, axes, angles = [], [], []
    named = dict(zip('xyz', np.eye(3)))
    for time, axis, angle in args.pulse:
        try:
            time, angle = float(time), float(angle)
            vector = named[axis] if axis in named else np.array([float(x) for x in axis.split(',')])
        except ValueError:
            parser.error('pulse times/angles must be numeric and axes must be x, y, z, or three numbers')
        if not np.isfinite(time) or time < 0 or not np.isfinite(angle):
            parser.error('pulse time must be finite/nonnegative and angle must be finite')
        if times and time < times[-1]:
            parser.error('pulse times must be nondecreasing; coincident pulses keep their order')
        if vector.shape != (3,) or not np.all(np.isfinite(vector)) or not np.linalg.norm(vector) > 0:
            parser.error('axis must be a finite nonzero three-component vector')
        times.append(time)
        angles.append(angle)
        axes.append(vector / np.linalg.norm(vector))
    # Exclusive creation prevents accidentally overwriting a sequence.
    with h5py.File(args.output, 'x') as file:
        group = file.create_group('pulses')
        group.create_dataset('times', data=np.array(times, dtype=float))
        group.create_dataset('axes', data=np.array(axes, dtype=float).reshape(-1, 3))
        group.create_dataset('angles', data=np.array(angles, dtype=float))
        group.attrs['rotation'] = 'positive +i angle axis.S; hbar=1; radians; measurement-grid pulses use Hamiltonian fields -angle/h*axis'


if __name__ == '__main__':
    main()
