#!/usr/bin/env python3
"""Plot Re m^z(t) for all Keldysh type-C 1M-sample datasets."""

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "Data" / "Keldysh_stype_C_hz0.5"
PLOTS_DIR = ROOT / "Plots"
SAMPLES_PER_CORE = 1_000_000


def load_dataset(path: Path) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5file:
        parameters = h5file["parameters"].attrs
        beta = float(parameters["beta"])
        samples_per_core = int(parameters["num_SamplesPerCore"])
        if samples_per_core != SAMPLES_PER_CORE:
            raise ValueError(f"Unexpected sample count in {path}: {samples_per_core}")

        delta_real_t = float(parameters["delta_real_t"])
        values = np.asarray(h5file["results/Re_magnetization"][:, 2], dtype=float)
        stds = np.asarray(
            h5file["runtimedata/Re_magnetization_sample_stds"][:, 2], dtype=float
        )

    times = delta_real_t * np.arange(values.size, dtype=float)
    return beta, times, values, stds


def main() -> None:
    plt.style.use(ROOT / "matplotlibrc")
    paths = sorted(DATA_DIR.glob(f"*samples_per_core={SAMPLES_PER_CORE}.hdf5"))
    if not paths:
        raise FileNotFoundError(f"No {SAMPLES_PER_CORE} samples-per-core files in {DATA_DIR}")

    datasets = sorted((load_dataset(path) for path in paths), key=lambda item: item[0])
    reference_times = datasets[0][1]
    if any(not np.array_equal(times, reference_times) for _, times, _, _ in datasets[1:]):
        raise ValueError("The selected beta datasets do not share a common real-time grid")

    figure, axis = plt.subplots(figsize=(6.5, 4.2))
    for beta, times, values, stds in datasets:
        (line,) = axis.plot(times, values, linewidth=1.5, label=rf"$\beta J_Q={beta:g}$")
        axis.fill_between(
            times,
            values - stds,
            values + stds,
            color=line.get_color(),
            alpha=0.20,
            linewidth=0.0,
        )

    axis.set_xlim(float(reference_times[0]), float(reference_times[-1]))
    axis.set_xlabel(r"$t J_Q$")
    axis.set_ylabel(r"$\mathrm{Re}\,m^z(t)$")
    axis.legend(loc="best", ncol=2)
    figure.tight_layout()

    PLOTS_DIR.mkdir(exist_ok=True)
    output = PLOTS_DIR / "Keldysh_stype_C_hz0.5_magnetization_1m"
    figure.savefig(f"{output}.png", dpi=300, bbox_inches="tight")
    figure.savefig(f"{output}.pdf", bbox_inches="tight")
    plt.close(figure)
    print(f"Saved {output}.png and {output}.pdf")


if __name__ == "__main__":
    main()
