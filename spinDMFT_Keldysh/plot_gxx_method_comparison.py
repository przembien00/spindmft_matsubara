#!/usr/bin/env python3
"""Compare the four available h_z=0 g^xx(t) data sources."""

from __future__ import annotations

import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data"
OUTPUT = ROOT / "Plots" / "gxx_method_comparison"
BETAS = (0.5, 1.5, 2.5)
TMAX = 15.0

KELDYSH_RE = re.compile(
    r"spinmodel=ISO__beta=(?P<beta>[0-9.]+)_samples_per_core=(?P<samples>[0-9]+)\.hdf5$"
)


def highest_keldysh_files() -> dict[float, Path]:
    selected: dict[float, tuple[int, Path]] = {}
    for path in (DATA / "Keldysh_stype_A").glob("*.hdf5"):
        match = KELDYSH_RE.fullmatch(path.name)
        if match is None:
            continue
        beta, samples = float(match["beta"]), int(match["samples"])
        if beta not in selected or samples > selected[beta][0]:
            selected[beta] = samples, path
    return {beta: path for beta, (_, path) in selected.items()}


def load_contour(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5file:
        if "real_times" in h5file["results"].attrs:
            times = h5file["results"].attrs["real_times"]
        else:
            parameters = h5file["parameters"].attrs
            times = np.linspace(0.0, parameters["Tmax"], parameters["num_RealTimePoints"])
        real = h5file["results/Re_correlation"][:, 0, 0]
        # Match the sign convention used by the analytic-continuation cxx data.
        imag = -h5file["results/Im_correlation"][:, 0, 0]
        real_err = h5file["runtimedata/Re_correlation_sample_stds"][:, 0, 0]
        imag_err = h5file["runtimedata/Im_correlation_sample_stds"][:, 0, 0]
    return times, real, imag, real_err, imag_err


def load_random(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5file:
        parameters = h5file["parameters"].attrs
        times = np.linspace(0.0, parameters["Tmax"], parameters["num_TimePoints"])
        real = h5file["results/Re_correlation"][0]
        imag = -h5file["results/Im_correlation"][0]
        real_err = h5file["results/Re_stddev"][0]
        imag_err = h5file["results/Im_stddev"][0]
    return times, real, imag, real_err, imag_err


def load_ana_cont(
    path: Path, beta: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5file:
        group = h5file[f"beta={beta:g}"]
        values = group["cxx"][...]
        lower = group["cxx_band_low"][...]
        upper = group["cxx_band_high"][...]
        return (
            group["t"][...], values.real, values.imag,
            lower.real, upper.real, lower.imag, upper.imag,
        )


def crop(times: np.ndarray, *values: np.ndarray) -> tuple[np.ndarray, ...]:
    mask = times <= TMAX + 1.0e-12
    return (times[mask], *(value[mask] for value in values))


def make_plot(
    real_data: dict[float, dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]],
    imag_data: dict[float, dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]],
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(10.0, 6.8), sharex="col", sharey="row")
    colors = {
        "Keldysh spinDMFT": "forestgreen",
        "random extrapolation": "crimson",
        "Unitary spinDMFT": "royalblue",
        "analytic continuation": "black",
    }
    styles = {
        "Keldysh spinDMFT": "-",
        "random extrapolation": "--",
        "Unitary spinDMFT": ":",
        "analytic continuation": "-.",
    }

    for row, (data, component) in enumerate(((real_data, r"\mathrm{Re}"), (imag_data, r"\mathrm{Im}"))):
        for axis, beta in zip(axes[row], BETAS):
            for method, (times, values, lower, upper) in data[beta].items():
                axis.fill_between(times, lower, upper, color=colors[method], alpha=0.18,
                                  linewidth=0, zorder=1)
                axis.plot(
                    times,
                    values,
                    color=colors[method],
                    linestyle=styles[method],
                    linewidth=1.7,
                    label=method,
                )
            if row == 0:
                axis.set_title(rf"$\beta J_Q={beta:g}$")
            axis.set_xlim(0.0, TMAX)
            if row == 1:
                axis.set_xlabel(r"$t J_Q$")
            if axis is axes[row, 0]:
                axis.set_ylabel(rf"${component}\,g^{{xx}}(t)$")
            axis.grid(True, alpha=0.35)

    axes[0, 0].legend(
        handles=[Line2D([], [], color=colors[name], linestyle=styles[name], label=name)
                 for name in colors],
        loc="best",
        fontsize="small",
    )
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.08)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        fig.savefig(OUTPUT / f"gxx_comparison_hz0.{extension}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    plt.style.use(ROOT / "matplotlibrc")
    keldysh_files = highest_keldysh_files()
    random_dir = DATA / "Benchmarks" / "Random_extrapolation"
    unitary_dir = DATA / "Benchmarks" / "Unitary_real_time"
    ana_path = DATA / "Benchmarks" / "spinDMFT_ana_cont" / "realtime_JL0_hz0_all_betas.hdf5"

    data: dict[float, dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]] = {}
    for beta in BETAS:
        k_times, k_real, k_imag, k_real_err, k_imag_err = crop(*load_contour(keldysh_files[beta]))
        r_times, r_real, r_imag, r_real_err, r_imag_err = crop(*load_random(
            random_dir / f"ISO__Random__N=inf__beta={beta:g}.hdf5"))
        u_times, u_real, u_imag, u_real_err, u_imag_err = crop(*load_contour(
            unitary_dir / f"spinmodel=ISO__beta={beta:g}.hdf5"))
        (a_times, a_real, a_imag, a_real_low, a_real_high,
         a_imag_low, a_imag_high) = crop(*load_ana_cont(ana_path, beta))
        data[beta] = {
            "Keldysh spinDMFT": (k_times, k_real, k_real - k_real_err, k_real + k_real_err),
            "random extrapolation": (r_times, r_real, r_real - r_real_err, r_real + r_real_err),
            "Unitary spinDMFT": (u_times, u_real, u_real - u_real_err, u_real + u_real_err),
            "analytic continuation": (a_times, a_real, np.minimum(a_real_low, a_real_high),
                                       np.maximum(a_real_low, a_real_high)),
        }
        data[beta + 1000] = {
            "Keldysh spinDMFT": (k_times, -k_imag, -k_imag - k_imag_err, -k_imag + k_imag_err),
            "random extrapolation": (r_times, -r_imag, -r_imag - r_imag_err, -r_imag + r_imag_err),
            "Unitary spinDMFT": (u_times, -u_imag, -u_imag - u_imag_err, -u_imag + u_imag_err),
            "analytic continuation": (a_times, -a_imag, -np.maximum(a_imag_low, a_imag_high),
                                       -np.minimum(a_imag_low, a_imag_high)),
        }

    make_plot(
        {beta: data[beta] for beta in BETAS},
        {beta: data[beta + 1000] for beta in BETAS},
    )


if __name__ == "__main__":
    main()
