#!/usr/bin/env python3
"""Fit and plot high-field correlations with Gaussian damping."""

from __future__ import annotations

import csv
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "Data" / "highB_highT"
OUTPUT_DIR = ROOT / "Plots" / "highB_highT_gaussian_fits"
FIELDS = (3, 7, 10)
INITIAL_CORRELATION = 0.25


def load_component(path: Path, component: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5file:
        results = h5file["results"]
        labels = list(results["correlation_direction_labels"].asstr()[...])
        if component not in labels:
            raise ValueError(f"{component!r} is not available in {path}")
        component_index = labels.index(component)
        times = np.asarray(results.attrs["real_times"], dtype=float)
        # The real-time correlation is stored at the tau=0 edge.
        values = np.asarray(results["Re_correlation"][:, component_index, 0], dtype=float)
        errors = np.asarray(
            h5file["runtimedata/Re_correlation_sample_stds"][:, component_index, 0],
            dtype=float,
        )
    return times, values, errors


def model_xx(times: np.ndarray, decay_time: float, frequency: float) -> np.ndarray:
    return np.exp(-(times / decay_time) ** 2) * np.cos(frequency * times)


def model_xy(times: np.ndarray, decay_time: float, frequency: float) -> np.ndarray:
    return np.exp(-(times / decay_time) ** 2) * np.sin(frequency * times)


def fit_component(
    times: np.ndarray, values: np.ndarray, errors: np.ndarray, component: str, field: int
) -> tuple[float, float]:
    # Normalize the measured correlations so the requested unit-amplitude forms
    # apply. Keep B signed: the stored xy convention is negative for positive h_z.
    normalized = values / INITIAL_CORRELATION
    normalized_errors = errors / INITIAL_CORRELATION
    positive_errors = normalized_errors[np.isfinite(normalized_errors) & (normalized_errors > 0)]
    if positive_errors.size == 0:
        raise ValueError("No positive finite correlation uncertainties are available")
    sigma = np.maximum(normalized_errors, np.median(positive_errors))

    model = model_xx if component == "xx" else model_xy
    frequency_guess = float(field if component == "xx" else -field)
    low_frequency, high_frequency = sorted((0.5 * frequency_guess, 1.5 * frequency_guess))
    parameters, _ = curve_fit(
        model,
        times,
        normalized,
        p0=(2.0, frequency_guess),
        sigma=sigma,
        absolute_sigma=True,
        bounds=((0.05, low_frequency), (100.0, high_frequency)),
        maxfev=100_000,
    )
    return float(parameters[0]), float(parameters[1])


def make_plot(component: str, datasets: dict[int, tuple[np.ndarray, ...]]) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8), sharey=True)
    for axis, field in zip(axes, FIELDS):
        times, values, errors, decay_time, frequency = datasets[field]
        data_color = "red"
        axis.fill_between(
            times,
            values - errors,
            values + errors,
            color=data_color,
            alpha=0.18,
            linewidth=0,
        )
        axis.plot(times, values, color=data_color, linewidth=1.35, label="data")

        dense_times = np.linspace(float(times[0]), float(times[-1]), 1600)
        fit_model = model_xx if component == "xx" else model_xy
        fit_values = INITIAL_CORRELATION * fit_model(dense_times, decay_time, frequency)
        axis.plot(dense_times, fit_values, color="blue", linestyle="--", linewidth=1.5,
                  label="fit")
        axis.text(
            0.04,
            0.96,
            f"$T={decay_time:.3g}$\n$B={frequency:.4g}$",
            transform=axis.transAxes,
            va="top",
            ha="left",
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.85,
                  "edgecolor": "0.75"},
        )
        axis.set_title(rf"$h_z={field}$")
        axis.set_xlim(float(times[0]), min(float(times[-1]), 10.0))
        axis.set_xlabel(r"$t J_Q$")
        axis.grid(True, alpha=0.35)

    axes[0].set_ylabel(rf"$\mathrm{{Re}}\,g^{{{component}}}(t)$")
    axes[0].legend(loc="lower left")
    fig.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output = OUTPUT_DIR / f"Re_g{component}_highB_highT"
    fig.savefig(f"{output}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{output}.pdf", bbox_inches="tight")
    plt.close(fig)
    return output


def main() -> None:
    plt.style.use(ROOT / "matplotlibrc")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fit_results: dict[str, dict[int, tuple[float, float]]] = {"xx": {}, "xy": {}}

    for component in ("xx", "xy"):
        datasets: dict[int, tuple[np.ndarray, ...]] = {}
        for field in FIELDS:
            path = DATA_DIR / f"spinmodel=ISO__beta=0.1__h=z_h_abs={field}X.hdf5"
            if not path.exists():
                raise FileNotFoundError(path)
            times, values, errors = load_component(path, component)
            if values[0] <= 0 and component == "xx":
                raise ValueError(f"Cannot normalize {component} in {path}: Re g^xx(0) <= 0")
            decay_time, frequency = fit_component(times, values, errors, component, field)
            fit_results[component][field] = decay_time, frequency
            datasets[field] = times, values, errors, decay_time, frequency
        print(make_plot(component, datasets))

    csv_path = OUTPUT_DIR / "fit_parameters.csv"
    with csv_path.open("w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(("component", "h_z", "T", "B", "fit_model", "source_file"))
        for component, field_results in fit_results.items():
            for field, (decay_time, frequency) in field_results.items():
                model = "exp(-(t/T)^2) cos(B t)" if component == "xx" else "exp(-(t/T)^2) sin(B t)"
                writer.writerow((component, field, f"{decay_time:.10g}",
                                 f"{frequency:.10g}", model,
                                 f"spinmodel=ISO__beta=0.1__h=z_h_abs={field}X.hdf5"))
    print(csv_path)


if __name__ == "__main__":
    main()
