#!/usr/bin/env python3
"""Compare highest-sample Keldysh real-time correlators with ana-cont data."""

from __future__ import annotations

import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data"
OUTPUT = ROOT / "Plots" / "Keldysh_ana_cont_comparison"
SAMPLE_RE = re.compile(r"beta=(?P<beta>[0-9.]+).*samples_per_core=(?P<samples>[0-9]+)\.hdf5$")


def highest_sample_files(directory: Path) -> dict[float, Path]:
    selected: dict[float, tuple[int, Path]] = {}
    for path in directory.glob("*.hdf5"):
        match = SAMPLE_RE.search(path.name)
        if match is None:
            continue
        beta, samples = float(match["beta"]), int(match["samples"])
        if beta not in selected or samples > selected[beta][0]:
            selected[beta] = (samples, path)
    return {beta: item[1] for beta, item in selected.items()}


def load_keldysh(path: Path, component: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5file:
        labels = list(h5file["results/correlation_direction_labels"].asstr()[...])
        if component == "zz" and "zz" not in labels and labels == ["xx"]:
            # At h_z=0, the ISO state is rotationally invariant.
            component = "xx"
        if component == "xy" and labels == ["xx"]:
            times = h5file["results"].attrs["real_times"]
            zeros = np.zeros_like(times, dtype=float)
            return times, zeros, zeros, zeros, zeros
        index = labels.index(component)
        times = h5file["results"].attrs["real_times"]
        return (
            times,
            h5file["results/Re_correlation"][:, index, 0],
            h5file["runtimedata/Re_correlation_sample_stds"][:, index, 0],
            h5file["results/Im_correlation"][:, index, 0],
            h5file["runtimedata/Im_correlation_sample_stds"][:, index, 0],
        )


def continuation_component(group: h5py.Group, component: str, zero_field: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if zero_field:
        if component == "xy":
            t = group["t"][...]
            return t, np.zeros_like(t, dtype=complex), np.zeros_like(t)
        central = group["cxx"][...]
        bands = group["xx_spectral"]
    elif component in {"xx", "xy"}:
        central = group[f"c{component}"][...]
        bands = group["transverse"]
    else:
        central = group["czz"][...]
        bands = group["longitudinal"]
    lower, upper = bands["band_low"][...], bands["band_high"][...]
    if lower.ndim == 2:
        transverse_index = {"xx": 0, "xy": 1}[component]
        lower, upper = lower[:, transverse_index], upper[:, transverse_index]
    error = np.maximum(np.abs(central - lower), np.abs(upper - central))
    return group["t"][...], central, error


def plot_case(keldysh_dir: Path, continuation_path: Path, case_stem: str, zero_field: bool) -> None:
    files = highest_sample_files(keldysh_dir)
    with h5py.File(continuation_path, "r") as continuation:
        betas = sorted(set(files).intersection(float(key.split("=")[1]) for key in continuation))
        for component in ("xx", "xy", "zz"):
            fig, axes = plt.subplots(2, 1, figsize=(8.0, 7.2), sharex=True, squeeze=False)
            axes = axes[:, 0]
            colors = plt.get_cmap("viridis")(np.linspace(0.10, 0.90, len(betas)))
            max_time = 0.0
            for beta, color in zip(betas, colors):
                t_k, re_k, re_err, im_k, im_err = load_keldysh(files[beta], component)
                group = continuation[f"beta={beta:g}"]
                t_a, ana, ana_err = continuation_component(group, component, zero_field)
                max_time = max(max_time, float(t_k[-1]), float(t_a[-1]))
                label = rf"$\beta={beta:g}$"
                for axis, k_value, k_error, a_value in (
                    (axes[0], re_k, re_err, ana.real),
                    (axes[1], -im_k, im_err, ana.imag),
                ):
                    axis.errorbar(t_a, a_value, yerr=ana_err.real if axis is axes[0] else ana_err.imag,
                                  color=color, linestyle="--", linewidth=1.2, errorevery=20,
                                  capsize=1.5, alpha=0.9, label=label, marker='^', markevery=20)
                    axis.errorbar(t_k, k_value, yerr=k_error, color=color, linestyle="-", linewidth=1.6,
                                  errorevery=8, capsize=1.5, alpha=0.95)
            axes[0].set_ylabel(rf"$\mathrm{{Re}}\,g^{{{component}}}(t)$")
            axes[1].set_ylabel(rf"$\mathrm{{Im}}\,g^{{{component}}}(t)$")
            axes[1].set_xlabel(r"$t J_Q$")
            axes[1].set_xlim(0.0, max_time)
            beta_legend = axes[0].legend(loc="best", ncol=2, title=r"$\beta$")
            axes[0].add_artist(beta_legend)
            axes[0].legend(handles=[
                Line2D([], [], color="black", linestyle="-", label="Keldysh"),
                Line2D([], [], color="black", linestyle="--", label="ana-cont"),
            ], loc="upper right")
            fig.tight_layout()
            for extension in ("png", "pdf"):
                fig.savefig(OUTPUT / f"{case_stem}_{component}.{extension}", dpi=300, bbox_inches="tight")
            plt.close(fig)


def main() -> None:
    plt.style.use(ROOT / "matplotlibrc")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plot_case(DATA / "Keldysh_stype_A", DATA / "spinDMFT_ana_cont/realtime_JL0_hz0_all_betas.hdf5", "hz0", True)
    plot_case(DATA / "Keldysh_stype_C_hz0.5", DATA / "spinDMFT_ana_cont/realtime_JL0_hz0.5_all_betas.hdf5", "hz0.5", False)


if __name__ == "__main__":
    main()
