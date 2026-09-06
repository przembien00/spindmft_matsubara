#!/usr/bin/env python3
"""Plot maximum real correlation and magnetization stds for a Keldysh data set."""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


FILENAME_RE = re.compile(r"__beta=(?P<beta>[0-9.]+).*_samples_per_core=(?P<samples>\d+)\.hdf5$")


def load_maximum_stds(data_dir: Path, dataset: str) -> dict[float, list[tuple[int, float]]]:
    data: dict[float, list[tuple[int, float]]] = defaultdict(list)
    for path in data_dir.glob("*.hdf5"):
        match = FILENAME_RE.search(path.name)
        if match is None:
            continue
        with h5py.File(path, "r") as h5file:
            samples_per_core = int(match["samples"])
            num_cores = int(h5file["parameters"].attrs["num_Cores"])
            data[float(match["beta"])].append(
                (samples_per_core * num_cores, float(np.max(h5file[dataset][...])))
            )
    return {beta: sorted(points) for beta, points in sorted(data.items())}


def plot_series(
    data: dict[float, list[tuple[int, float]]], ylabel: str, output: Path
) -> None:
    fig, axis = plt.subplots(figsize=(6.5, 4.2))
    all_stds = []
    for beta, points in data.items():
        samples, stds = np.asarray(points, dtype=float).T
        all_stds.extend(stds)
        axis.plot(samples, stds, marker="o", label=rf"$\beta={beta:g}$")

    axis.set_xscale("log")
    if any(std > 0.0 for std in all_stds):
        axis.set_yscale("log")
    else:
        axis.axhline(0.0, color="black", linewidth=0.8, zorder=0)
        axis.set_ylim(-0.05, 0.05)
    axis.set_xlabel("Total samples")
    axis.set_ylabel(ylabel)
    axis.legend(ncol=2, loc="best")

    output.parent.mkdir(exist_ok=True)
    fig.savefig(Path(f"{output}.pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    args = parser.parse_args()
    data_dir = args.data_dir.resolve()
    plt.style.use(Path(__file__).resolve().parent / "matplotlibrc")
    plots_dir = Path(__file__).resolve().parent / "Plots"
    output_prefix = data_dir.name
    plot_series(
        load_maximum_stds(data_dir, "runtimedata/Re_correlation_sample_stds"),
        r"$\max\, \sigma[\mathrm{Re}\, g^{ab}_{ij}]$",
        plots_dir / f"{output_prefix}_largest_real_correlation_std_vs_total_samples",
    )
    plot_series(
        load_maximum_stds(data_dir, "runtimedata/Re_magnetization_sample_stds"),
        r"$\max\, \sigma[\mathrm{Re}\, m^a]$",
        plots_dir / f"{output_prefix}_largest_real_magnetization_std_vs_total_samples",
    )


if __name__ == "__main__":
    main()
