#!/usr/bin/env python3
"""Compare h_z=0.5 xx, zz, and xy correlations across three betas."""

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
COMPONENT_INDEX = {"xx": 0, "xy": 1, "zz": 3}
KELDYSH_RE = re.compile(
    r"spinmodel=ISO__beta=(?P<beta>[0-9.]+)__h=z_h_abs=0.5_samples_per_core=(?P<samples>[0-9]+)\.hdf5$"
)


def highest_keldysh_files() -> dict[float, Path]:
    selected: dict[float, tuple[int, Path]] = {}
    for path in (DATA / "Keldysh_stype_C_hz0.5").glob("*.hdf5"):
        match = KELDYSH_RE.fullmatch(path.name)
        if match is None:
            continue
        beta, samples = float(match["beta"]), int(match["samples"])
        if beta not in selected or samples > selected[beta][0]:
            selected[beta] = samples, path
    return {beta: path for beta, (_, path) in selected.items()}


def crop(times: np.ndarray, *values: np.ndarray) -> tuple[np.ndarray, ...]:
    mask = times <= TMAX + 1.0e-12
    return (times[mask], *(value[mask] for value in values))


def load_contour(path: Path, component: str) -> tuple[np.ndarray, ...]:
    with h5py.File(path, "r") as h5file:
        labels = [label.decode() for label in h5file["results/correlation_direction_labels"][...]]
        index = labels.index(component)
        times = h5file["results"].attrs.get("real_times")
        if times is None:
            parameters = h5file["parameters"].attrs
            times = np.linspace(0.0, parameters["Tmax"], parameters["num_RealTimePoints"])
        real = h5file["results/Re_correlation"][:, index, 0]
        imag = h5file["results/Im_correlation"][:, index, 0]
        real_err = h5file["runtimedata/Re_correlation_sample_stds"][:, index, 0]
        imag_err = h5file["runtimedata/Im_correlation_sample_stds"][:, index, 0]
    return times, real, imag, real_err, imag_err


def load_random(path: Path, component: str) -> tuple[np.ndarray, ...]:
    with h5py.File(path, "r") as h5file:
        parameters = h5file["parameters"].attrs
        times = np.linspace(0.0, parameters["Tmax"], parameters["num_TimePoints"])
        index = COMPONENT_INDEX[component]
        real = h5file["results/Re_correlation"][index]
        imag = h5file["results/Im_correlation"][index]
        real_err = h5file["results/Re_stddev"][index]
        imag_err = h5file["results/Im_stddev"][index]
    return times, real, imag, real_err, imag_err


def load_ana_cont(path: Path, beta: float, component: str) -> tuple[np.ndarray, ...]:
    with h5py.File(path, "r") as h5file:
        group = h5file[f"beta={beta:g}"]
        if component == "zz":
            values = group["czz"][...]
            bands = group["longitudinal"]
            lower = bands["band_low"][...]
            upper = bands["band_high"][...]
        else:
            values = group[f"c{component}"][...]
            # The transverse band_* arrays are not the envelopes of cxx/cxy
            # in the finite-field continuation output.  Use the matching
            # continuation ensemble directly instead.
            ensemble = group[f"c{component}_ensemble"][...]
            lower = np.min(ensemble.real, axis=0) + 1j * np.min(ensemble.imag, axis=0)
            upper = np.max(ensemble.real, axis=0) + 1j * np.max(ensemble.imag, axis=0)
        return group["t"][...], values.real, values.imag, lower.real, upper.real, lower.imag, upper.imag


def make_plot(component: str, data: dict[float, dict[str, tuple[np.ndarray, ...]]]) -> None:
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
    for row, (part, label) in enumerate((("real", r"\mathrm{Re}"), ("imag", r"\mathrm{Im}"))):
        for axis, beta in zip(axes[row], BETAS):
            for method, (times, values, lower, upper) in data[beta][part].items():
                axis.fill_between(times, lower, upper, color=colors[method], alpha=0.18,
                                  linewidth=0, zorder=1)
                axis.plot(times, values, color=colors[method], linestyle=styles[method],
                          linewidth=1.7, label=method, zorder=2)
            if row == 0:
                axis.set_title(rf"$\beta J_Q={beta:g}$")
            axis.set_xlim(0.0, TMAX)
            axis.grid(True, alpha=0.35)
            if row == 1:
                axis.set_xlabel(r"$t J_Q$")
            if axis is axes[row, 0]:
                axis.set_ylabel(rf"${label}\,g^{{{component}}}(t)$")
    axes[0, 0].legend(
        handles=[Line2D([], [], color=colors[name], linestyle=styles[name], label=name)
                 for name in colors], loc="best", fontsize="small",
    )
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.08)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        fig.savefig(OUTPUT / f"g{component}_comparison_hz0p5.{extension}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    plt.style.use(ROOT / "matplotlibrc")
    keldysh_files = highest_keldysh_files()
    random_dir = DATA / "Benchmarks" / "Random_extrapolation"
    unitary_dir = DATA / "Benchmarks" / "Unitary_real_time"
    ana_path = DATA / "Benchmarks" / "spinDMFT_ana_cont" / "realtime_JL0_hz0.5_all_betas.hdf5"

    for component in ("xx", "zz", "xy"):
        data: dict[float, dict[str, dict[str, tuple[np.ndarray, ...]]]] = {}
        for beta in BETAS:
            k = crop(*load_contour(keldysh_files[beta], component))
            r = crop(*load_random(
                random_dir / f"ISO__Random__N=inf__beta={beta:g}__h_z=0.5.hdf5", component))
            a = crop(*load_ana_cont(ana_path, beta, component))
            random_real = -r[1] if component == "xy" else r[1]
            random_imag = -r[2] if component == "xy" else r[2]
            data[beta] = {"real": {}, "imag": {}}
            data[beta]["real"]["Keldysh spinDMFT"] = (k[0], k[1], k[1] - k[3], k[1] + k[3])
            data[beta]["real"]["random extrapolation"] = (r[0], random_real, random_real - r[3], random_real + r[3])
            data[beta]["real"]["analytic continuation"] = (a[0], a[1], np.minimum(a[3], a[4]), np.maximum(a[3], a[4]))
            data[beta]["imag"]["Keldysh spinDMFT"] = (k[0], k[2], k[2] - k[4], k[2] + k[4])
            data[beta]["imag"]["random extrapolation"] = (r[0], random_imag, random_imag - r[4], random_imag + r[4])
            data[beta]["imag"]["analytic continuation"] = (a[0], -a[2], -np.maximum(a[5], a[6]), -np.minimum(a[5], a[6]))
            unitary_path = unitary_dir / f"spinmodel=ISO__beta={beta:g}__h=z_h_abs=0.5.hdf5"
            if unitary_path.exists():
                u = crop(*load_contour(unitary_path, component))
                data[beta]["real"]["Unitary spinDMFT"] = (u[0], u[1], u[1] - u[3], u[1] + u[3])
                data[beta]["imag"]["Unitary spinDMFT"] = (u[0], u[2], u[2] - u[4], u[2] + u[4])
        make_plot(component, data)


if __name__ == "__main__":
    main()
