#!/usr/bin/env python3
"""Plot the JL=0 zz correlations at beta=0.5 and beta=2.5."""

from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data" / "Keldysh_stype_C_hz0.5"
LOW_T_DATA = ROOT / "Data" / "Keldysh_stype_C_hz0.5_discretization_separate"
OUTPUT = ROOT / "Plots" / "Keldysh_stype_C_hz0.5_gzz_temperature_comparison"
SAMPLES_PER_CORE = 500_000
LOW_T_PATH = LOW_T_DATA / (
    "spinmodel=ISO__beta=2.5__h=z_h_abs=0.5_scan=real_time__"
    "numImagTimeSteps=80__numRealTimeSteps=400.hdf5"
)


def load_correlation(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5file:
        results = h5file["results"]
        labels = results["correlation_direction_labels"].asstr()[...].tolist()
        component = labels.index("zz")
        times = results.attrs["real_times"]
        values = results["Re_correlation"][:, component, 0]
    return times, values


def main() -> None:
    plt.style.use(ROOT / "matplotlibrc")
    fig, ax = plt.subplots(figsize=(6.4, 4.1))

    curves = (
        (
            DATA / f"spinmodel=ISO__beta=0.5__h=z_h_abs=0.5_samples_per_core={SAMPLES_PER_CORE}.hdf5",
            "red",
            "high T",
            6.3,
            8,
        ),
        (LOW_T_PATH, "blue", "low T", 6.3, 8),
    )
    for path, color, label, label_time, label_offset in curves:
        times, values = load_correlation(path)
        ax.plot(times, values, color=color, linewidth=1.8)
        index = int(np.argmin(np.abs(times - label_time)))
        ax.annotate(
            label,
            xy=(times[index], values[index]),
            xytext=(0, label_offset),
            textcoords="offset points",
            color=color,
            ha="center",
            va="bottom",
        )

    ax.set_xlim(0.0, 15.0)
    ax.set_xlabel(r"$t J_Q$")
    ax.set_ylabel(r"$\mathrm{Re}\,g^{zz}(t)$")
    fig.tight_layout()

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        # Build the filename directly: with_suffix() treats the decimal in
        # "hz0.5" as a suffix and would drop the rest of the output stem.
        fig.savefig(Path(f"{OUTPUT}.{extension}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
