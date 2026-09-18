#!/usr/bin/env python3
"""Plot the finite-field unitary real-time longitudinal correlation."""

from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "Data" / "Benchmarks" / "Unitary_real_time"
OUTPUT_DIR = ROOT / "Plots"


def load_gzz(path: Path) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """Return beta and the tau=0 real-time gzz components from one file."""
    with h5py.File(path, "r") as h5file:
        parameters = h5file["parameters"].attrs
        beta = float(parameters["beta"])
        labels = list(h5file["results/correlation_direction_labels"].asstr()[...])
        zz_index = labels.index("zz")
        times = np.linspace(0.0, parameters["Tmax"], parameters["num_RealTimePoints"])
        real = h5file["results/Re_correlation"][:, zz_index, 0]
        imag = h5file["results/Im_correlation"][:, zz_index, 0]
    return beta, times, real, imag


def main() -> None:
    plt.style.use(ROOT / "matplotlibrc")
    curves = sorted(
        (load_gzz(path) for path in DATA_DIR.glob("spinmodel=ISO__beta=*__h=z_h_abs=0.5.hdf5")),
        key=lambda curve: curve[0],
    )
    if not curves:
        raise FileNotFoundError(f"No finite-field unitary files found in {DATA_DIR}")

    fig, axes = plt.subplots(2, 1, figsize=(6.4, 6.6), sharex=True)
    for beta, times, real, imag in curves:
        label = rf"$\beta J_Q={beta:g}$"
        axes[0].plot(times, real, label=label)
        axes[1].plot(times, imag, label=label)

    axes[0].set_ylabel(r"$\mathrm{Re}\,g^{zz}(t)$")
    axes[1].set_ylabel(r"$\mathrm{Im}\,g^{zz}(t)$")
    axes[1].set_xlabel(r"$t J_Q$")
    axes[0].legend(loc="best", ncol=2)
    for axis in axes:
        axis.set_xlim(curves[0][1][0], curves[0][1][-1])

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.08)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        fig.savefig(OUTPUT_DIR / f"unitary_gzz_all_betas_hz0p5.{extension}", dpi=300,
                    bbox_inches="tight")


if __name__ == "__main__":
    main()
