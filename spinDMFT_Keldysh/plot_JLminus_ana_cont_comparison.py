#!/usr/bin/env python3
"""Compare negative-JL Keldysh and two analytic-continuation datasets."""

from __future__ import annotations

import re
import colorsys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import numpy as np
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data"
OUTPUT = ROOT / "Plots" / "JLminus_ana_cont_comparison"
BETAS = (0.5, 1.5, 2.5)
TMAX = 15.0
MAGNETIZATION_LIGHTNESS = 0.58
MAGNETIZATION_SATURATION = 0.95
COMPONENT_INDEX = {"xx": 0, "xy": 1, "zz": 3}
KELDYSH_RE = re.compile(
    r"spinmodel=ISO__JL=-2\.4495__beta=(?P<beta>[0-9.]+)__h=z_h_abs=0\.5_"
    r"samples_per_core=(?P<samples>[0-9]+)\.hdf5$"
)


def highest_keldysh_files() -> dict[float, Path]:
    selected: dict[float, tuple[int, Path]] = {}
    directory = DATA / "Keldysh_stype_C_hz0.5_JLm2p4495"
    for path in directory.glob("*.hdf5"):
        match = KELDYSH_RE.fullmatch(path.name)
        if match is None:
            continue
        beta, samples = float(match["beta"]), int(match["samples"])
        if beta not in BETAS:
            continue
        if beta not in selected or samples > selected[beta][0]:
            selected[beta] = samples, path
    missing = set(BETAS) - set(selected)
    if missing:
        raise FileNotFoundError(f"Missing Keldysh files for beta={sorted(missing)}")
    return {beta: path for beta, (_, path) in selected.items()}


def crop(times: np.ndarray, *arrays: np.ndarray) -> tuple[np.ndarray, ...]:
    mask = times <= TMAX + 1.0e-12
    return (times[mask], *(array[mask] for array in arrays))


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


def load_keldysh(path: Path, component: str) -> tuple[np.ndarray, ...]:
    with h5py.File(path, "r") as h5file:
        labels = [label.decode() if isinstance(label, bytes) else str(label)
                  for label in h5file["results/correlation_direction_labels"][...]]
        index = labels.index(component)
        results = h5file["results"]
        times = results.attrs.get("real_times")
        if times is None:
            parameters = h5file["parameters"].attrs
            times = np.linspace(0.0, parameters["Tmax"], parameters["num_RealTimePoints"])
        values = results["Re_correlation"][:, index, 0] + 1j * results["Im_correlation"][:, index, 0]
        errors = h5file["runtimedata"]
        re_err = errors["Re_correlation_sample_stds"][:, index, 0]
        im_err = errors["Im_correlation_sample_stds"][:, index, 0]
    return times, values.real, values.imag, re_err, im_err


def load_continuation(path: Path, beta: float, component: str) -> tuple[np.ndarray, ...]:
    with h5py.File(path, "r") as h5file:
        group = h5file[f"beta={beta:g}"]
        key = f"c{component}"
        values = group[f"{key}_re"][...] + 1j * group[f"{key}_im"][...]
        if component == "zz":
            ensemble_group = group["longitudinal"]
            ensemble_re = ensemble_group["real_time_ensemble_re"][...]
            ensemble_im = ensemble_group["real_time_ensemble_im"][...]
        else:
            ensemble_re = group[f"{key}_ensemble_re"][...]
            ensemble_im = group[f"{key}_ensemble_im"][...]
        lower_re, upper_re = np.min(ensemble_re, axis=0), np.max(ensemble_re, axis=0)
        lower_im, upper_im = np.min(ensemble_im, axis=0), np.max(ensemble_im, axis=0)
        times = group["t"][...]
    return times, values.real, values.imag, lower_re, upper_re, lower_im, upper_im


def load_keldysh_magnetization(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5file:
        times = h5file["results"].attrs.get("real_times")
        if times is None:
            parameters = h5file["parameters"].attrs
            times = np.linspace(0.0, parameters["Tmax"], parameters["num_RealTimePoints"])
        magnetization = h5file["results/Re_magnetization"][:, 2]
    return times, magnetization


def load_continuation_magnetization(path: Path, beta: float, *, qmc: bool) -> float:
    with h5py.File(path, "r") as h5file:
        group = h5file[f"beta={beta:g}"]
        if qmc:
            # QMC's equal-time transverse commutator gives m_z=2 Im gxy(0).
            return float(2.0 * group["cxy_im"][0])
        return float(group.attrs["magnetization_z"])


def plot_component(component: str, keldysh_files: dict[float, Path],
                   spin_cont_path: Path, qmc_cont_path: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(10.0, 7.0), sharex="col", sharey="row")
    colors = {
        "Keldysh spinDMFT": "forestgreen",
        "spinDMFT ana-cont": "black",
        "QMC ana-cont": "darkorange",
    }
    styles = {"Keldysh spinDMFT": "-", "spinDMFT ana-cont": "-.", "QMC ana-cont": "--"}

    for column, beta in enumerate(BETAS):
        kt, kre, kim, kre_err, kim_err = crop(*load_keldysh(keldysh_files[beta], component))
        series: dict[str, dict[str, tuple[np.ndarray, ...]]] = {
            "Keldysh spinDMFT": {
                "real": (kt, kre, kre - kre_err, kre + kre_err),
                "imag": (kt, kim, kim - kim_err, kim + kim_err),
            }
        }
        magnetization_products: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        if component == "zz":
            kt_m, km = load_keldysh_magnetization(keldysh_files[beta])
            kt_m, km = crop(kt_m, km)
            magnetization_products["Keldysh spinDMFT"] = (kt_m, km * km[0])
        for label, path in (("spinDMFT ana-cont", spin_cont_path), ("QMC ana-cont", qmc_cont_path)):
            at, are, aim, rlo, rhi, ilo, ihi = crop(*load_continuation(path, beta, component))
            # The continuation files use the opposite imaginary-part sign from
            # the stored Keldysh correlation, as in the existing h_z=0.5 plots.
            series[label] = {
                "real": (at, are, np.minimum(rlo, rhi), np.maximum(rlo, rhi)),
                "imag": (at, -aim, -np.maximum(ilo, ihi), -np.minimum(ilo, ihi)),
            }
            if component == "zz":
                m_z = load_continuation_magnetization(
                    path, beta, qmc=(label == "QMC ana-cont")
                )
                magnetization_products[label] = (at, np.full_like(at, m_z * m_z))

        for row, part in enumerate(("real", "imag")):
            axis = axes[row, column]
            for method, parts in series.items():
                times, values, lower, upper = parts[part]
                axis.fill_between(times, lower, upper, color=colors[method], alpha=0.18,
                                  linewidth=0, zorder=1)
                axis.plot(times, values, color=colors[method], linestyle=styles[method],
                          linewidth=1.7, label=method, zorder=2)
            if component == "zz" and row == 0:
                for method, (times, product) in magnetization_products.items():
                    axis.plot(times, product, color=vivid_color(colors[method]), linestyle=":",
                              linewidth=2.1, zorder=3)
            if row == 0:
                axis.set_title(rf"$\beta J_Q={beta:g}$")
            if row == 1:
                axis.set_xlabel(r"$t J_Q$")
            axis.set_xlim(0.0, TMAX)
            axis.grid(True, alpha=0.35)
            if column == 0:
                axis.set_ylabel(rf"${'\\mathrm{Re}' if part == 'real' else '\\mathrm{Im}'}\,g^{{{component}}}(t)$")

    if component == "zz":
        legend_handles = [
            Line2D([], [], color=colors[name], linestyle=styles[name], label=name)
            for name in colors
        ] + [
            Line2D([], [], color="black", linestyle=":",
                   label=rf"$\langle S^z(t)\rangle\langle S^z(0)\rangle$")
        ]
        axes[0, 0].legend(legend_handles,
                          [handle.get_label() for handle in legend_handles],
                          loc="upper left", fontsize="small")
        fig.tight_layout()
    else:
        axes[0, 0].legend(
            handles=[Line2D([], [], color=colors[name], linestyle=styles[name], label=name)
                     for name in colors],
            loc="best", fontsize="small",
        )
        fig.tight_layout()
    fig.subplots_adjust(hspace=0.08)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        fig.savefig(OUTPUT / f"g{component}_comparison_JLminus_hz0p5.{extension}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    plt.style.use(ROOT / "matplotlibrc")
    keldysh_files = highest_keldysh_files()
    spin_cont_path = DATA / "Benchmarks" / "spinDMFT_ana_cont" / "realtime_JLminus_hz0.5_all_betas.hdf5"
    qmc_cont_path = DATA / "Benchmarks" / "QMC_ana_cont" / "realtime_Jminus_hz0.5_N32768_all_betas.hdf5"
    for path in (spin_cont_path, qmc_cont_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    for component in ("xx", "xy", "zz"):
        plot_component(component, keldysh_files, spin_cont_path, qmc_cont_path)


if __name__ == "__main__":
    main()
