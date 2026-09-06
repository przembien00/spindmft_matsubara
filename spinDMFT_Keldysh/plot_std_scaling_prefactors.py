#!/usr/bin/env python3
"""Fit max real-component uncertainties to C / sqrt(N) and plot C(beta)."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


KELDYSH_DIR = Path(__file__).resolve().parent
DATA_DIR = KELDYSH_DIR / "Data"
PLOTS_DIR = KELDYSH_DIR / "Plots"
FILENAME_RE = re.compile(r"__beta=(?P<beta>[0-9.]+).*_samples_per_core=(?P<samples>\d+)\.hdf5$")
CASES = (
    ("Keldysh_stype_A", r"A, $h_z=0$", "o"),
    ("Keldysh_stype_A_antithetic", r"A antithetic, $h_z=0$", "s"),
    ("Keldysh_stype_C_hz0.5", r"C, $h_z=0.5$", "^"),
    ("Keldysh_stype_C_hz0.5_antithetic", r"C antithetic, $h_z=0.5$", "D"),
)


def load_maximum_stds(directory: Path, dataset: str) -> dict[float, list[tuple[int, float]]]:
    """Load (total trajectory count, maximum real std) grouped by beta."""
    data: dict[float, list[tuple[int, float]]] = defaultdict(list)
    for path in directory.glob("*.hdf5"):
        match = FILENAME_RE.search(path.name)
        if match is None:
            continue
        with h5py.File(path, "r") as h5file:
            samples_per_core = int(match["samples"])
            num_cores = int(h5file["parameters"].attrs["num_Cores"])
            total_samples = samples_per_core * num_cores
            data[float(match["beta"])].append(
                (total_samples, float(np.max(h5file[dataset][...])))
            )
    return {beta: sorted(points) for beta, points in sorted(data.items())}


def fit_prefactors(data: dict[float, list[tuple[int, float]]]) -> dict[float, float]:
    """Constrained least-squares fit of std(N) = C / sqrt(N) for each beta."""
    prefactors = {}
    for beta, points in data.items():
        samples, stds = np.asarray(points, dtype=float).T
        inverse_sqrt_samples = 1.0 / np.sqrt(samples)
        prefactors[beta] = float(
            np.linalg.lstsq(inverse_sqrt_samples[:, None], stds, rcond=None)[0][0]
        )
    return prefactors


def plot_prefactors(dataset: str, ylabel: str, output_name: str, *, skip_zero: bool) -> None:
    fig, axis = plt.subplots(figsize=(6.5, 4.2))
    for directory_name, label, marker in CASES:
        prefactors = fit_prefactors(load_maximum_stds(DATA_DIR / directory_name, dataset))
        betas = np.asarray(sorted(prefactors), dtype=float)
        values = np.asarray([prefactors[beta] for beta in betas])
        if skip_zero and np.all(values == 0.0):
            continue
        axis.plot(betas, values, marker=marker, label=label)

    axis.set_xlabel(r"$\beta$")
    axis.set_ylabel(ylabel)
    axis.legend(loc="best")
    fig.savefig(PLOTS_DIR / output_name, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    plt.style.use(KELDYSH_DIR / "matplotlibrc")
    plot_prefactors(
        "runtimedata/Re_correlation_sample_stds",
        r"$C_g$ in $\max\,\sigma[\mathrm{Re}\,g^{ab}_{ij}] = C_g/\sqrt{N}$",
        "largest_real_correlation_std_scaling_prefactor_vs_beta.pdf",
        skip_zero=False,
    )
    plot_prefactors(
        "runtimedata/Re_magnetization_sample_stds",
        r"$C_m$ in $\max\,\sigma[\mathrm{Re}\,m^a] = C_m/\sqrt{N}$",
        "largest_real_magnetization_std_scaling_prefactor_vs_beta.pdf",
        skip_zero=True,
    )


if __name__ == "__main__":
    main()
