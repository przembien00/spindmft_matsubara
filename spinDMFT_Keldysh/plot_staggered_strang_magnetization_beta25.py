#!/usr/bin/env python3
"""Plot beta=2.5 longitudinal magnetization for the two hz=0.5 staged cases."""

from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data" / "Staggered Strang staged"
OUTPUT = ROOT / "Plots" / "Staggered_Strang_ana_cont_comparison"
CASES = (
    ("JL0_hz0p5", r"$J_L=0$",
     DATA / "staged__spinmodel=ISO__beta=2.5__h=z_h_abs=0.5__prop=staggered-strang.hdf5",
     "forestgreen"),
    ("JLminus_hz0p5", r"$J_L=-2.4495$",
     DATA / "staged__spinmodel=ISO__JL=-2.4495__beta=2.5__h=z_h_abs=0.5__prop=staggered-strang.hdf5",
     "darkorange"),
)


def load_magnetization(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5file:
        times = h5file["results"].attrs["real_times"]
        magnetization_z = h5file["results/Re_magnetization"][:, 2]
        std_z = h5file["runtimedata/Re_magnetization_sample_stds"][:, 2]
    return times, magnetization_z, std_z


def main() -> None:
    plt.style.use(ROOT / "matplotlibrc")
    for case, label, path, color in CASES:
        fig, axis = plt.subplots(figsize=(8.0, 4.8))
        times, magnetization, std = load_magnetization(path)
        initial_magnetization = magnetization[0]
        axis.fill_between(times, initial_magnetization - std, initial_magnetization + std,
                          color=color, alpha=0.22, linewidth=0, label="_nolegend_")
        axis.plot(times, magnetization, color=color, linewidth=1.8,
                  label=label, zorder=3)
        envelope_min = float(min(np.min(magnetization), np.min(initial_magnetization - std)))
        envelope_max = float(max(np.max(magnetization), np.max(initial_magnetization + std)))
        margin = max(0.12 * (envelope_max - envelope_min), 1.0e-5)
        axis.set_ylim(envelope_min - margin, envelope_max + margin)
        axis.set_xlim(0.0, 15.0)
        axis.set_xlabel(r"$t J_Q$")
        axis.set_ylabel(r"$\langle S^z(t)\rangle$")
        axis.grid(True, alpha=0.3)
        fig.tight_layout()
        case_output = OUTPUT / case
        case_output.mkdir(parents=True, exist_ok=True)
        fig.savefig(case_output / "magnetization_beta2p5.pdf", bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
