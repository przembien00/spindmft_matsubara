#!/usr/bin/env python3
"""Plot the real-time z-magnetization range versus real-time step size."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DISCRETIZATION_DIR = ROOT / "Data" / "Keldysh_stype_C_hz0.5_discretization"
BASELINE_DIR = ROOT / "Data" / "Keldysh_stype_C_hz0.5"
PLOTS_DIR = ROOT / "Plots"
SAMPLES_PER_CORE = 100_000
ALTERNATE_BETA = 2.5
ALTERNATE_SAMPLES_PER_CORE = 1_000_000
BASELINE_STEPS = 200


def read_point(path: Path) -> tuple[float, float, int, int, float, float]:
    """Return beta, step size, grid/sample counts, range, and propagated error."""
    with h5py.File(path, "r") as h5file:
        parameters = h5file["parameters"].attrs
        beta = float(parameters["beta"])
        delta_real_t = float(parameters["delta_real_t"])
        num_real_steps = int(parameters["num_RealTimeSteps"])
        samples_per_core = int(parameters["num_SamplesPerCore"])
        magnetization_z = h5file["results/Re_magnetization"][:, 2]
        magnetization_z_stds = h5file["runtimedata/Re_magnetization_sample_stds"][:, 2]
    minimum_index = int(np.argmin(magnetization_z))
    maximum_index = int(np.argmax(magnetization_z))
    magnetization_range = float(magnetization_z[maximum_index] - magnetization_z[minimum_index])
    magnetization_range_std = float(
        np.hypot(magnetization_z_stds[maximum_index], magnetization_z_stds[minimum_index])
    )
    return beta, delta_real_t, num_real_steps, samples_per_core, magnetization_range, magnetization_range_std


def load_points() -> dict[float, list[tuple[float, int, int, float, float]]]:
    points: dict[float, list[tuple[float, int, int, float, float]]] = defaultdict(list)

    for path in DISCRETIZATION_DIR.glob("*.hdf5"):
        beta, delta_real_t, num_real_steps, samples_per_core, magnetization_range, magnetization_range_std = read_point(path)
        if samples_per_core != SAMPLES_PER_CORE:
            raise ValueError(f"Unexpected sample count in {path}")
        points[beta].append(
            (delta_real_t, num_real_steps, samples_per_core, magnetization_range, magnetization_range_std)
        )

    for path in BASELINE_DIR.glob("*.hdf5"):
        beta, delta_real_t, num_real_steps, samples_per_core, magnetization_range, magnetization_range_std = read_point(path)
        expected_samples = (
            ALTERNATE_SAMPLES_PER_CORE if beta == ALTERNATE_BETA else SAMPLES_PER_CORE
        )
        if num_real_steps == BASELINE_STEPS and samples_per_core == expected_samples:
            points[beta].append(
                (delta_real_t, num_real_steps, samples_per_core, magnetization_range, magnetization_range_std)
            )

    expected_betas = (0.2, 0.5, 1.0, 1.5, 2.0, 2.5)
    for beta in expected_betas:
        if len(points[beta]) < 2:
            raise ValueError(
                f"Expected at least two discretizations for beta={beta:g}; "
                f"found {len(points[beta])} points."
            )
    return {beta: sorted(values) for beta, values in sorted(points.items())}


def main() -> None:
    plt.style.use(ROOT / "matplotlibrc")
    points = load_points()

    fig, axis = plt.subplots(figsize=(6.5, 4.2))
    for beta, values in points.items():
        delta_real_t, _, _, magnetization_range, magnetization_range_std = np.asarray(values, dtype=float).T
        axis.errorbar(
            delta_real_t,
            magnetization_range,
            yerr=magnetization_range_std,
            fmt="o-",
            capsize=3,
            label=rf"$\beta J_\mathrm{{Q}}={beta:g}$",
        )

    all_steps = [step for values in points.values() for step, _, _, _, _ in values]
    all_lower_bounds = [
        magnetization_range - magnetization_range_std
        for values in points.values()
        for _, _, _, magnetization_range, magnetization_range_std in values
    ]
    axis.set_xlim(min(all_steps), max(all_steps))
    axis.set_ylim(bottom=min(all_lower_bounds))
    axis.set_xlabel(r"Real-time step size $\Delta t_\mathrm{R}$")
    axis.set_ylabel(r"$\max_t\,\mathrm{Re}\,m^z(t)-\min_t\,\mathrm{Re}\,m^z(t)$")
    axis.legend(ncol=2, loc="best")

    PLOTS_DIR.mkdir(exist_ok=True)
    output = PLOTS_DIR / "hz0p5_magnetization_z_range_vs_real_step_size"
    fig.savefig(f"{output}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{output}.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
