#!/usr/bin/env python3
"""Compare h_z=0.5 xx, zz, and xy correlations across three betas."""

from __future__ import annotations

import re
import colorsys
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import numpy as np
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data"
OUTPUT = ROOT / "Plots" / "gxx_method_comparison"
BETAS = (0.5, 1.5, 2.5)
TMAX = 15.0
COMPONENT_INDEX = {"xx": 0, "xy": 1, "zz": 3}
MAGNETIZATION_LIGHTNESS = 0.58
MAGNETIZATION_SATURATION = 0.95
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


def vivid_color(color: str) -> tuple[float, float, float]:
    """Brighten mean-product colors while preserving their hue and saturation."""
    rgb = to_rgb(color)
    hue, lightness, saturation = colorsys.rgb_to_hls(*rgb)
    if saturation == 0.0:
        return (0.40, 0.40, 0.40)
    return colorsys.hls_to_rgb(
        hue,
        max(lightness, MAGNETIZATION_LIGHTNESS),
        max(saturation, MAGNETIZATION_SATURATION),
    )


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


def load_contour_magnetization(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5file:
        results = h5file["results"]
        times = results.attrs.get("real_times")
        if times is None:
            parameters = h5file["parameters"].attrs
            times = np.linspace(0.0, parameters["Tmax"], parameters["num_RealTimePoints"])
        magnetization = results["Re_magnetization"][:, 2]
    return times, magnetization


def load_initial_magnetization_from_gxy(path: Path) -> float:
    """Use the equal-time identity m_z = 2 Im g^xy(0)."""
    with h5py.File(path, "r") as h5file:
        return float(2.0 * h5file["results/Im_correlation"][COMPONENT_INDEX["xy"], 0])


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
        name = f"c{component}"
        # The current finite-field continuation schema stores the central
        # values and ensembles as separate real/imaginary datasets.
        if name in group:
            values = group[name][...]
        else:
            values = group[f"{name}_re"][...] + 1j * group[f"{name}_im"][...]

        if component == "zz":
            bands = group["longitudinal"]
            if "band_low" in bands:
                lower, upper = bands["band_low"][...], bands["band_high"][...]
            else:
                lower = bands["band_low_re"][...] + 1j * bands["band_low_im"][...]
                upper = bands["band_high_re"][...] + 1j * bands["band_high_im"][...]
        elif f"{name}_ensemble" in group:
            ensemble = group[f"{name}_ensemble"][...]
            lower = np.min(ensemble.real, axis=0) + 1j * np.min(ensemble.imag, axis=0)
            upper = np.max(ensemble.real, axis=0) + 1j * np.max(ensemble.imag, axis=0)
        else:
            ensemble = group[f"{name}_ensemble_re"][...] + 1j * group[f"{name}_ensemble_im"][...]
            lower = np.min(ensemble.real, axis=0) + 1j * np.min(ensemble.imag, axis=0)
            upper = np.max(ensemble.real, axis=0) + 1j * np.max(ensemble.imag, axis=0)
        return group["t"][...], values.real, values.imag, lower.real, upper.real, lower.imag, upper.imag


def make_plot(
    component: str,
    data: dict[float, dict[str, tuple[np.ndarray, ...]]],
    magnetization_products: dict[float, dict[str, tuple[np.ndarray, np.ndarray]]],
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(10.0, 6.8), sharex="col", sharey="row")
    colors = {
        "Keldysh spinDMFT": "forestgreen",
        "Chebyshev": "crimson",
        "analytic continuation": "black",
    }
    styles = {
        "Keldysh spinDMFT": "-",
        "Chebyshev": "--",
        "analytic continuation": "-.",
    }
    for row, (part, label) in enumerate((("real", r"\mathrm{Re}"), ("imag", r"\mathrm{Im}"))):
        for column, (axis, beta) in enumerate(zip(axes[row], BETAS)):
            for method, (times, values, lower, upper) in data[beta][part].items():
                axis.fill_between(times, lower, upper, color=colors[method], alpha=0.18,
                                  linewidth=0, zorder=1)
                axis.plot(times, values, color=colors[method], linestyle=styles[method],
                          linewidth=1.7, label=method, zorder=2)
            if component == "zz" and row == 0:
                for method, (times, product) in magnetization_products[beta].items():
                    axis.plot(
                        times,
                        product,
                        color=vivid_color(colors[method]),
                        linestyle=":",
                        linewidth=2.1,
                        label="_nolegend_",
                        zorder=3,
                    )
            if row == 0:
                axis.set_title(rf"$\beta J_Q={beta:g}$")
            axis.set_xlim(0.0, TMAX)
            axis.grid(True, alpha=0.35)
            if row == 1:
                axis.set_xlabel(r"$t J_Q$")
            if axis is axes[row, 0]:
                axis.set_ylabel(rf"${label}\,g^{{{component}}}(t)$")
    handles = [Line2D([], [], color=colors[name], linestyle=styles[name], label=name)
               for name in colors]
    if component == "zz":
        handles.append(Line2D(
            [], [], color="black", linestyle=":",
            label=r"$\langle S^z(t)\rangle\langle S^z(0)\rangle$",
        ))
    axes[0, 0].legend(handles=handles, loc="best", fontsize="small")
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
    ana_path = DATA / "Benchmarks" / "spinDMFT_ana_cont" / "realtime_JL0_hz0.5_all_betas.hdf5"

    for component in ("xx", "zz", "xy"):
        data: dict[float, dict[str, dict[str, tuple[np.ndarray, ...]]]] = {}
        magnetization_products: dict[float, dict[str, tuple[np.ndarray, np.ndarray]]] = {}
        for beta in BETAS:
            k = crop(*load_contour(keldysh_files[beta], component))
            r = crop(*load_random(
                random_dir / f"ISO__Random__N=inf__beta={beta:g}__h_z=0.5.hdf5", component))
            a = crop(*load_ana_cont(ana_path, beta, component))
            random_real = -r[1] if component == "xy" else r[1]
            random_imag = -r[2] if component == "xy" else r[2]
            data[beta] = {"real": {}, "imag": {}}
            data[beta]["real"]["Keldysh spinDMFT"] = (k[0], k[1], k[1] - k[3], k[1] + k[3])
            data[beta]["real"]["Chebyshev"] = (r[0], random_real, random_real - r[3], random_real + r[3])
            data[beta]["real"]["analytic continuation"] = (a[0], a[1], np.minimum(a[3], a[4]), np.maximum(a[3], a[4]))
            data[beta]["imag"]["Keldysh spinDMFT"] = (k[0], k[2], k[2] - k[4], k[2] + k[4])
            data[beta]["imag"]["Chebyshev"] = (r[0], random_imag, random_imag - r[4], random_imag + r[4])
            data[beta]["imag"]["analytic continuation"] = (a[0], -a[2], -np.maximum(a[5], a[6]), -np.minimum(a[5], a[6]))

            if component == "zz":
                mag_times, mag_k = crop(*load_contour_magnetization(keldysh_files[beta]))
                magnetization_products[beta] = {
                    "Keldysh spinDMFT": (mag_times, mag_k * mag_k[0]),
                }
                random_path = random_dir / f"ISO__Random__N=inf__beta={beta:g}__h_z=0.5.hdf5"
                m_cheb = load_initial_magnetization_from_gxy(random_path)
                cheb_times, = crop(r[0])
                magnetization_products[beta]["Chebyshev"] = (
                    cheb_times, np.full_like(cheb_times, m_cheb * m_cheb)
                )
                with h5py.File(ana_path, "r") as h5file:
                    spin_ac_m = float(h5file[f"beta={beta:g}"].attrs["magnetization_z"])
                magnetization_products[beta]["analytic continuation"] = (
                    a[0], np.full_like(a[0], spin_ac_m * spin_ac_m)
                )
        make_plot(component, data, magnetization_products)


if __name__ == "__main__":
    main()
