#!/usr/bin/env python3
"""Compare init-therm and Keldysh real-time greater/lesser correlations."""

from __future__ import annotations

import re
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parent
INIT_DATA = ROOT.parent / "spinDMFT_init_therm" / "Data"
KELDYSH_DATA = ROOT / "Data"
BETAS = (0.5, 1.5, 2.5)
COMPONENTS = ("xx", "xy", "yx", "zz")
TMAX = 15.0

KELDYSH_PATTERNS = {
    "stype_A": re.compile(
        r"spinmodel=ISO__beta=(?P<beta>[0-9.]+)_samples_per_core=(?P<samples>[0-9]+)\.hdf5$"
    ),
    "stype_C_hz0p5": re.compile(
        r"spinmodel=ISO__beta=(?P<beta>[0-9.]+)__h=z_h_abs=0\.5_samples_per_core=(?P<samples>[0-9]+)\.hdf5$"
    ),
}
KELDYSH_DIRECTORY_NAMES = {
    "stype_A": "Keldysh_stype_A",
    "stype_C_hz0p5": "Keldysh_stype_C_hz0.5",
}


def select_keldysh_files(dataset: str) -> dict[float, Path]:
    selected: dict[float, tuple[int, Path]] = {}
    directory = KELDYSH_DATA / KELDYSH_DIRECTORY_NAMES[dataset]
    pattern = KELDYSH_PATTERNS[dataset]
    for path in directory.glob("*.hdf5"):
        match = pattern.fullmatch(path.name)
        if match is None:
            continue
        beta, samples = float(match["beta"]), int(match["samples"])
        if beta in BETAS and (beta not in selected or samples > selected[beta][0]):
            selected[beta] = (samples, path)
    missing = set(BETAS) - set(selected)
    if missing:
        raise FileNotFoundError(f"Missing {dataset} Keldysh files for beta={sorted(missing)}")
    return {beta: value[1] for beta, value in selected.items()}


def init_file(beta: float, field: bool) -> Path:
    field_tag = "__h=z_h_abs=0.5" if field else ""
    path = INIT_DATA / f"spinmodel=ISO__beta={beta:g}{field_tag}__init_term.hdf5"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def decode_labels(values: np.ndarray) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def load_init(path: Path, component: str) -> tuple[np.ndarray, dict[str, tuple[np.ndarray, np.ndarray]]]:
    with h5py.File(path, "r") as data:
        directions = data["parameters/correlation_directions"][...].astype(int)
        labels = ["".join(("xyz"[a], "xyz"[b])) for a, b in directions]
        index = labels.index(component)
        times = data["results/real_times"][...]
        values = {}
        endpoint = "greater"
        real = data[f"results/Re_{endpoint}"][index]
        imag = data[f"results/Im_{endpoint}"][index]
        real_err = data[f"runtimedata/Re_{endpoint}_sample_stds"][index]
        imag_err = data[f"runtimedata/Im_{endpoint}_sample_stds"][index]
        values[endpoint] = (real + 1j * imag, real_err + 1j * imag_err)
    return times, values


def load_keldysh(path: Path, component: str) -> tuple[np.ndarray, dict[str, tuple[np.ndarray, np.ndarray]]]:
    with h5py.File(path, "r") as data:
        labels = decode_labels(data["results/correlation_direction_labels"][...])
        index = labels.index(component)
        times = data["results"].attrs["real_times"]
        values = {}
        # The edge-grid contract is tau=beta -> greater.
        for endpoint, edge in (("greater", -1),):
            real = data["results/Re_correlation"][:, index, edge]
            imag = data["results/Im_correlation"][:, index, edge]
            real_err = data["runtimedata/Re_correlation_sample_stds"][:, index, edge]
            imag_err = data["runtimedata/Im_correlation_sample_stds"][:, index, edge]
            values[endpoint] = (real + 1j * imag, real_err + 1j * imag_err)
    return times, values


def load_init_magnetization_z(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as data:
        times = data["results/real_times"][...]
        values = data["results/Re_magnetization"][:, 0]
    return times, values


def load_keldysh_magnetization_z(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as data:
        times = data["results"].attrs["real_times"]
        values = data["results/Re_magnetization"][:, 2]
    return times, values


def plot_component(dataset: str, component: str, files: dict[float, Path], field: bool) -> None:
    output = ROOT / "Plots" / f"init_therm_vs_Keldysh_{dataset}"
    fig, axes = plt.subplots(2, 3, figsize=(10.6, 6.8), sharex="col", sharey="row")
    colors = {"init-therm": "royalblue", "Keldysh": "darkorange"}

    for col, beta in enumerate(BETAS):
        init_times, init_values = load_init(init_file(beta, field), component)
        keldysh_times, keldysh_values = load_keldysh(files[beta], component)
        if field and component == "zz":
            init_m_times, init_m = load_init_magnetization_z(init_file(beta, field))
            keldysh_m_times, keldysh_m = load_keldysh_magnetization_z(files[beta])
        for row, part in enumerate(("real", "imag")):
            for source, times, values in (
                ("init-therm", init_times, init_values),
                ("Keldysh", keldysh_times, keldysh_values),
            ):
                for endpoint in ("greater",):
                    correlation, error = values[endpoint]
                    central = correlation.real if part == "real" else correlation.imag
                    sigma = error.real if part == "real" else error.imag
                    axes[row, col].fill_between(
                        times, central - sigma, central + sigma,
                        color=colors[source], alpha=0.11, linewidth=0,
                    )
                    axes[row, col].plot(
                        times, central, color=colors[source], linewidth=1.65,
                    )
            if field and component == "zz" and row == 0:
                for source, times, magnetization in (
                    ("init-therm", init_m_times, init_m),
                    ("Keldysh", keldysh_m_times, keldysh_m),
                ):
                    axes[row, col].plot(
                        times, magnetization * magnetization[0],
                        color=colors[source], linestyle=":", linewidth=1.8,
                    )
            axes[row, col].set_xlim(0.0, TMAX)
            axes[row, col].grid(True, alpha=0.35)
            if row == 0:
                axes[row, col].set_title(rf"$\beta J_Q={beta:g}$")
            if row == 1:
                axes[row, col].set_xlabel(r"$t J_Q$")
            if col == 0:
                axes[row, col].set_ylabel(
                    rf"$\mathrm{{{'Re' if part == 'real' else 'Im'}}}\,g^{{{component}}}(t)$"
                )

    if field and component == "zz":
        handles = [
            Line2D([], [], color=colors[source], linestyle=linestyle, linewidth=1.8,
                   label=label)
            for source, linestyle, label in (
                ("init-therm", "-", "init-therm greater"),
                ("Keldysh", "-", "Keldysh greater"),
                ("init-therm", ":", r"init-therm $m_z(t)m_z(0)$"),
                ("Keldysh", ":", r"Keldysh $m_z(t)m_z(0)$"),
            )
        ]
        fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
                   bbox_to_anchor=(0.5, -0.005))
        fig.tight_layout(rect=(0, 0.055, 1, 1))
    else:
        handles = [
            Line2D([], [], color=colors[source], linewidth=1.8, label=f"{source} greater")
            for source in ("init-therm", "Keldysh")
        ]
        fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
                   bbox_to_anchor=(0.5, -0.005))
        fig.tight_layout(rect=(0, 0.055, 1, 1))
    output.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        fig.savefig(output / f"g{component}_comparison.{extension}", dpi=300,
                    bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    plt.style.use(ROOT / "matplotlibrc")
    configurations = (("stype_A", False, ("xx",)), ("stype_C_hz0p5", True, COMPONENTS))
    for dataset, field, components in configurations:
        files = select_keldysh_files(dataset)
        print(f"{dataset} selected Keldysh files:")
        for beta in BETAS:
            print(f"  beta={beta:g}: {files[beta].name}")
        for component in components:
            plot_component(dataset, component, files, field)


if __name__ == "__main__":
    main()
