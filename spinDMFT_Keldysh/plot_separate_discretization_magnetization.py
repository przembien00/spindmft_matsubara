#!/usr/bin/env python3
"""Plot the real-time z-magnetization range for separate Keldysh grid scans."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = ROOT / "Data" / "Keldysh_stype_C_hz0.5_discretization_separate"
PLOTS_DIR = ROOT / "Plots"


def read_point(path: Path) -> tuple[str, float, float, float, int, int, float, float]:
    """Return scan type, beta, both grid spacings, grid counts, and magnetization range."""
    with h5py.File(path, "r") as h5file:
        parameters = h5file["parameters"].attrs
        extension = parameters["(original) filename_extension"].decode()
        scan_type = extension.split("scan=", 1)[1].split("__", 1)[0]
        beta = float(parameters["beta"])
        delta_imag_t = float(parameters["delta_imag_t"])
        delta_real_t = float(parameters["delta_real_t"])
        num_imag_steps = int(parameters["num_ImagTimeSteps"])
        num_real_steps = int(parameters["num_RealTimeSteps"])
        magnetization = h5file["results/Re_magnetization"][:, 2]
        stddev = h5file["runtimedata/Re_magnetization_sample_stds"][:, 2]

    low = int(np.argmin(magnetization))
    high = int(np.argmax(magnetization))
    span = float(magnetization[high] - magnetization[low])
    span_std = float(np.hypot(stddev[high], stddev[low]))
    return scan_type, beta, delta_imag_t, delta_real_t, num_imag_steps, num_real_steps, span, span_std


def load_points() -> dict[str, dict[float, list[tuple[float, float, float]]]]:
    points: dict[str, dict[float, list[tuple[float, float, float]]]] = {
        "imaginary_time": defaultdict(list),
        "real_time": defaultdict(list),
    }
    for path in sorted(DATA_DIR.glob("*.hdf5")):
        scan_type, beta, delta_imag_t, delta_real_t, num_imag_steps, num_real_steps, span, span_std = read_point(path)
        if scan_type not in points:
            raise ValueError(f"Unknown scan type {scan_type!r} in {path}")
        if scan_type == "imaginary_time":
            points["imaginary_time"][beta].append((delta_imag_t, span, span_std))
        else:
            points["real_time"][beta].append((delta_real_t, span, span_std))
            # The 200-step real-time run supplies the matching N_tau=80 point
            # for the imaginary-time scan while keeping the real-time grid fixed.
            if num_real_steps == 200 and num_imag_steps == 80:
                points["imaginary_time"][beta].append((delta_imag_t, span, span_std))

    expected_betas = (0.2, 0.5, 1.0, 1.5, 2.0, 2.5)
    for scan_type, beta_points in points.items():
        expected_count = 3
        for beta in expected_betas:
            if len(beta_points[beta]) != expected_count:
                raise ValueError(
                    f"Expected {expected_count} {scan_type} discretizations at beta={beta:g}; "
                    f"found {len(beta_points[beta])}."
                )
            beta_points[beta].sort()
    return points


def plot_scan(
    scan_type: str,
    points: dict[float, list[tuple[float, float, float]]],
    output_suffix: str,
) -> None:
    fig, axis = plt.subplots(figsize=(6.5, 4.2))
    for beta, values in sorted(points.items()):
        x, span, span_std = np.asarray(values, dtype=float).T
        axis.errorbar(
            x,
            span,
            yerr=span_std,
            fmt="o-",
            capsize=3,
            label=rf"$\beta J_\mathrm{{Q}}={beta:g}$",
        )

    all_x = [row[0] for values in points.values() for row in values]
    all_lows = [row[1] - row[2] for values in points.values() for row in values]
    axis.set_xlim(min(all_x), max(all_x))
    axis.set_ylim(bottom=min(all_lows))
    if scan_type == "imaginary_time":
        axis.set_xlabel(r"Imaginary-time step size $\Delta\tau$")
        output_stem = f"hz0p5_magnetization_z_range_vs_imaginary_step_size{output_suffix}"
    else:
        axis.set_xlabel(r"Real-time step size $\Delta t_\mathrm{R}$")
        output_stem = f"hz0p5_magnetization_z_range_vs_real_step_size_separate_scan{output_suffix}"
    axis.set_ylabel(r"$\max_t\,\mathrm{Re}\,m^z(t)-\min_t\,\mathrm{Re}\,m^z(t)$")
    axis.legend(ncol=2, loc="best")

    output = PLOTS_DIR / output_stem
    fig.savefig(f"{output}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{output}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-suffix", default="")
    args = parser.parse_args()

    plt.style.use(ROOT / "matplotlibrc")
    global DATA_DIR
    DATA_DIR = args.data_dir
    points = load_points()
    PLOTS_DIR.mkdir(exist_ok=True)
    for scan_type, beta_points in points.items():
        plot_scan(scan_type, beta_points, args.output_suffix)


if __name__ == "__main__":
    main()
