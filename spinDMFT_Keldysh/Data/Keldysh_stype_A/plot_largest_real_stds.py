#!/usr/bin/env python3
"""Plot the largest real-component uncertainty versus samples per MPI core."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


DATA_DIR = Path(__file__).resolve().parent
PLOTS_DIR = DATA_DIR / "Plots"
FILENAME_RE = re.compile(
    r"spinmodel=ISO__beta=(?P<beta>[0-9.]+)_samples_per_core=(?P<samples>\d+)\.hdf5"
)


def load_maximum_stds(dataset: str) -> dict[float, list[tuple[int, float]]]:
    """Return max(dataset) for every beta and samples-per-core input file."""
    data: dict[float, list[tuple[int, float]]] = defaultdict(list)
    for path in DATA_DIR.glob("*.hdf5"):
        match = FILENAME_RE.fullmatch(path.name)
        if match is None:
            continue
        with h5py.File(path, "r") as h5file:
            maximum_std = float(np.max(h5file[dataset][...]))
        data[float(match["beta"])].append((int(match["samples"]), maximum_std))
    return {beta: sorted(points) for beta, points in sorted(data.items())}


def plot_series(
    data: dict[float, list[tuple[int, float]]], ylabel: str, stem: str, *, log_y: bool
) -> None:
    fig, axis = plt.subplots(figsize=(6.5, 4.2))
    for beta, points in data.items():
        samples, stds = np.asarray(points, dtype=float).T
        axis.plot(samples, stds, marker="o", label=rf"$\beta={beta:g}$")

    axis.set_xscale("log")
    if log_y:
        axis.set_yscale("log")
    else:
        # All supplied real-magnetization standard-deviation arrays are zero.
        # A small symmetric range makes that fact visible without perturbing data.
        axis.axhline(0.0, color="black", linewidth=0.8, zorder=0)
        axis.set_ylim(-0.05, 0.05)
    axis.set_xlabel("Samples per MPI core")
    axis.set_ylabel(ylabel)
    axis.legend(ncol=2, loc="best")

    PLOTS_DIR.mkdir(exist_ok=True)
    for extension in ("png", "pdf"):
        fig.savefig(PLOTS_DIR / f"{stem}.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    plt.style.use(DATA_DIR.parents[1] / "matplotlibrc")
    plot_series(
        load_maximum_stds("runtimedata/Re_correlation_sample_stds"),
        r"$\max\, \sigma[\mathrm{Re}\, g^{ab}_{ij}]$",
        "largest_real_correlation_std_vs_samples_per_core",
        log_y=True,
    )
    plot_series(
        load_maximum_stds("runtimedata/Re_magnetization_sample_stds"),
        r"$\max\, \sigma[\mathrm{Re}\, m^a]$",
        "largest_real_magnetization_std_vs_samples_per_core",
        log_y=False,
    )


if __name__ == "__main__":
    main()
