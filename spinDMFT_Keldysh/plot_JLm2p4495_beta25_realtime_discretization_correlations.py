#!/usr/bin/env python3
"""Plot nonzero beta=2.5 correlations across real-time discretizations."""

from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "Data" / "Keldysh_stype_C_hz0.5_JLm2p4495_discretization_separate"
OUTPUT_DIR = ROOT / "Plots" / "JLm2p4495_beta2p5_realtime_discretizations"
BETWEEN_COMPONENTS = ("xx", "xy", "zz")
EXPECTED_STEPS = (100, 200, 400)
COLORS = {100: "royalblue", 200: "darkorange", 400: "forestgreen"}


def load_data() -> dict[int, dict[str, tuple[np.ndarray, ...]]]:
    """Read the tau=0 correlation slices for the three real-time grids."""
    data: dict[int, dict[str, tuple[np.ndarray, ...]]] = {}
    for path in sorted(DATA_DIR.glob("*beta=2.5*scan=real_time*.hdf5")):
        with h5py.File(path, "r") as h5file:
            parameters = h5file["parameters"].attrs
            if float(parameters["JL"]) >= 0:
                continue
            steps = int(parameters["num_RealTimeSteps"])
            labels = [label.decode() for label in h5file["results/correlation_direction_labels"][...]]
            real_times = h5file["results"].attrs["real_times"]
            if len(real_times) != int(parameters["num_RealTimePoints"]):
                raise ValueError(f"Real-time grid size mismatch in {path}")

            for component in BETWEEN_COMPONENTS:
                index = labels.index(component)
                values_re = h5file["results/Re_correlation"][:, index, 0]
                values_im = h5file["results/Im_correlation"][:, index, 0]
                errors_re = h5file["runtimedata/Re_correlation_sample_stds"][:, index, 0]
                errors_im = h5file["runtimedata/Im_correlation_sample_stds"][:, index, 0]
                data.setdefault(steps, {})[component] = (
                    np.asarray(real_times),
                    np.asarray(values_re),
                    np.asarray(errors_re),
                    np.asarray(values_im),
                    np.asarray(errors_im),
                )

    if tuple(sorted(data)) != EXPECTED_STEPS:
        raise ValueError(f"Expected real-time step counts {EXPECTED_STEPS}; found {tuple(sorted(data))}")
    for steps, components in data.items():
        missing = set(BETWEEN_COMPONENTS) - components.keys()
        if missing:
            raise ValueError(f"Missing components {sorted(missing)} for num_RealTimeSteps={steps}")
    return data


def main() -> None:
    plt.style.use(ROOT / "matplotlibrc")
    data = load_data()
    fig, axes = plt.subplots(3, 2, figsize=(9.0, 8.0), sharex="col")

    for row, component in enumerate(BETWEEN_COMPONENTS):
        for col, (part, part_label, title) in enumerate(
            (("real", r"\mathrm{Re}", "Real part"), ("imag", r"\mathrm{Im}", "Imaginary part"))
        ):
            axis = axes[row, col]
            for steps in EXPECTED_STEPS:
                times, values_re, errors_re, values_im, errors_im = data[steps][component]
                values = values_re if part == "real" else values_im
                errors = errors_re if part == "real" else errors_im
                axis.plot(times, values, color=COLORS[steps], linewidth=1.6)
                axis.fill_between(
                    times,
                    values - errors,
                    values + errors,
                    color=COLORS[steps],
                    alpha=0.14,
                    linewidth=0,
                )
            axis.set_ylabel(rf"${part_label}\,g^{{{component}}}(t,\tau=0)$")
            axis.grid(True, alpha=0.35)
            if row == 0:
                axis.set_title(title)
            if row == len(BETWEEN_COMPONENTS) - 1:
                axis.set_xlabel(r"Real time $t J_\mathrm{Q}$")

    handles = [
        Line2D([], [], color=COLORS[steps], linewidth=1.8, label=rf"$N_{{t}}={steps}$")
        for steps in EXPECTED_STEPS
    ]
    axes[0, 0].legend(handles=handles, loc="best", ncol=1)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.12, wspace=0.18)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        fig.savefig(OUTPUT_DIR / f"gxx_gxy_gzz_beta2p5_realtime_discretizations.{extension}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
