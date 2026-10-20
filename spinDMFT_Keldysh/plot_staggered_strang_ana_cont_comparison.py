#!/usr/bin/env python3
"""Compare staged staggered-Strang results with spinDMFT and QMC continuations."""

from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data"
STAGED_DIR = DATA / "Staggered Strang staged"
OUTPUT = ROOT / "Plots" / "Staggered_Strang_ana_cont_comparison"
PLOT_BETAS = (0.5, 1.5, 2.5)
RANDOM_COMPONENT_INDEX = {"xx": 0, "xy": 1, "zz": 3}
TMAX = 15.0


def staged_files(jl: float, h_abs: float, betas: tuple[float, ...]) -> dict[float, Path]:
    """Select staged data at the requested longitudinal coupling and field."""
    selected: dict[float, Path] = {}
    for path in sorted(STAGED_DIR.glob("*.hdf5")):
        with h5py.File(path, "r") as h5file:
            params = h5file["parameters"].attrs
            beta = float(params["beta"])
            if not np.isclose(float(params.get("JL", 0.0)), jl):
                continue
            if not np.isclose(float(params.get("h_abs", 0.0)), h_abs):
                continue
            if beta not in betas:
                continue
            if beta in selected:
                raise ValueError(f"Multiple matching staged files for beta={beta:g}")
            selected[beta] = path
    missing = set(betas) - set(selected)
    if missing:
        raise FileNotFoundError(f"Missing staged files for beta={sorted(missing)}")
    return selected


def crop(times: np.ndarray, *arrays: np.ndarray) -> tuple[np.ndarray, ...]:
    mask = times <= TMAX + 1.0e-12
    return (times[mask], *(array[mask] for array in arrays))


def load_staged(path: Path, component: str) -> tuple[np.ndarray, ...]:
    with h5py.File(path, "r") as h5file:
        results = h5file["results"]
        labels = [label.decode() if isinstance(label, bytes) else str(label)
                  for label in results["correlation_direction_labels"][...]]
        index = labels.index(component)
        times = results.attrs["real_times"]
        values = results["Re_correlation"][:, index, 0] + 1j * results["Im_correlation"][:, index, 0]
        re_err = h5file["runtimedata/Re_correlation_sample_stds"][:, index, 0]
        im_err = h5file["runtimedata/Im_correlation_sample_stds"][:, index, 0]
    return times, values.real, values.imag, re_err, im_err


def load_continuation(path: Path, beta: float, component: str) -> tuple[np.ndarray, ...]:
    with h5py.File(path, "r") as h5file:
        group = h5file[f"beta={beta:g}"]
        key = f"c{component}"
        if f"{key}_re" in group:
            values = group[f"{key}_re"][...] + 1j * group[f"{key}_im"][...]
            if component == "zz":
                ensemble = group["longitudinal"]
                ensemble_re = ensemble["real_time_ensemble_re"][...]
                ensemble_im = ensemble["real_time_ensemble_im"][...]
            else:
                ensemble_re = group[f"{key}_ensemble_re"][...]
                ensemble_im = group[f"{key}_ensemble_im"][...]
            lower_re, upper_re = np.min(ensemble_re, axis=0), np.max(ensemble_re, axis=0)
            lower_im, upper_im = np.min(ensemble_im, axis=0), np.max(ensemble_im, axis=0)
        else:
            values = group[key][...]
            lower, upper = group[f"{key}_band_low"][...], group[f"{key}_band_high"][...]
            lower_re, upper_re = np.minimum(lower.real, upper.real), np.maximum(lower.real, upper.real)
            lower_im, upper_im = np.minimum(lower.imag, upper.imag), np.maximum(lower.imag, upper.imag)
        times = group["t"][...]
    return times, values.real, values.imag, lower_re, upper_re, lower_im, upper_im


def load_magnetization(path: Path, beta: float, *, qmc: bool) -> float:
    with h5py.File(path, "r") as h5file:
        group = h5file[f"beta={beta:g}"]
        if qmc:
            return float(2.0 * group["cxy_im"][0])
        return float(group.attrs["magnetization_z"])


def load_staged_magnetization(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as h5file:
        results = h5file["results"]
        return results.attrs["real_times"], results["Re_magnetization"][:, 2]


def load_random(path: Path, component: str) -> tuple[np.ndarray, ...]:
    with h5py.File(path, "r") as h5file:
        parameters = h5file["parameters"].attrs
        times = np.linspace(0.0, parameters["Tmax"], parameters["num_TimePoints"])
        index = RANDOM_COMPONENT_INDEX[component]
        results = h5file["results"]
        return (
            times,
            results["Re_correlation"][index],
            results["Im_correlation"][index],
            results["Re_stddev"][index],
            results["Im_stddev"][index],
        )


def random_magnetization(path: Path) -> float:
    """Get m_z from the equal-time identity m_z=2 Im g^xy(0)."""
    with h5py.File(path, "r") as h5file:
        return float(2.0 * h5file["results/Im_correlation"][RANDOM_COMPONENT_INDEX["xy"], 0])


def plot_component(component: str, files: dict[float, Path], betas: tuple[float, ...],
                   spin_cont_path: Path, qmc_cont_path: Path | None,
                   random_h_abs: float | None, case_stem: str) -> None:
    fig, axes = plt.subplots(2, len(betas), figsize=(4.0 * len(betas), 7.0),
                             sharex="col", sharey="row", squeeze=False)
    colors = {
        "Staggered-Strang spinDMFT": "forestgreen",
        "Random-coupling extrapolation": "crimson",
        "spinDMFT ana-cont": "black",
        "QMC ana-cont": "darkorange",
    }
    styles = {
        "Staggered-Strang spinDMFT": "-",
        "Random-coupling extrapolation": "--",
        "spinDMFT ana-cont": "-.",
        "QMC ana-cont": "--",
    }

    methods = ["Staggered-Strang spinDMFT", "spinDMFT ana-cont"]
    if random_h_abs is not None:
        methods.insert(1, "Random-coupling extrapolation")
    if qmc_cont_path is not None:
        methods.append("QMC ana-cont")
    for column, beta in enumerate(betas):
        st, sre, sim, sre_err, sim_err = crop(*load_staged(files[beta], component))
        series: dict[str, dict[str, tuple[np.ndarray, ...]]] = {
            "Staggered-Strang spinDMFT": {
                "real": (st, sre, sre - sre_err, sre + sre_err),
                "imag": (st, sim, sim - sim_err, sim + sim_err),
            }
        }
        magnetization_products: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        if component == "zz":
            mt, magnetization = load_staged_magnetization(files[beta])
            mt, magnetization = crop(mt, magnetization)
            magnetization_products["Staggered-Strang spinDMFT"] = (mt, magnetization * magnetization[0])

        if random_h_abs is not None:
            suffix = "" if random_h_abs == 0.0 else f"__h_z={random_h_abs:g}"
            random_path = DATA / "Benchmarks" / "Random_extrapolation" / (
                f"ISO__Random__N=inf__beta={beta:g}{suffix}.hdf5"
            )
            rt, rr, ri, rr_err, ri_err = crop(*load_random(random_path, component))
            if component == "xy":
                rr, ri = -rr, -ri
            series["Random-coupling extrapolation"] = {
                "real": (rt, rr, rr - rr_err, rr + rr_err),
                "imag": (rt, ri, ri - ri_err, ri + ri_err),
            }
            if component == "zz":
                random_m = random_magnetization(random_path)
                magnetization_products["Random-coupling extrapolation"] = (
                    rt, np.full_like(rt, random_m * random_m)
                )

        continuation_paths = [("spinDMFT ana-cont", spin_cont_path)]
        if qmc_cont_path is not None:
            continuation_paths.append(("QMC ana-cont", qmc_cont_path))
        for label, path in continuation_paths:
            at, are, aim, rlo, rhi, ilo, ihi = crop(*load_continuation(path, beta, component))
            # Continuation files store the opposite imaginary-part sign to the real-time solver.
            series[label] = {
                "real": (at, are, np.minimum(rlo, rhi), np.maximum(rlo, rhi)),
                "imag": (at, -aim, -np.maximum(ilo, ihi), -np.minimum(ilo, ihi)),
            }
            if component == "zz":
                m_z = load_magnetization(path, beta, qmc=(label == "QMC ana-cont"))
                magnetization_products[label] = (at, np.full_like(at, m_z * m_z))

        for row, part in enumerate(("real", "imag")):
            axis = axes[row, column]
            for method, parts in series.items():
                times, values, lower, upper = parts[part]
                axis.fill_between(times, lower, upper, color=colors[method], alpha=0.18,
                                  linewidth=0, zorder=1)
                axis.plot(times, values, color=colors[method], linestyle=styles[method],
                          linewidth=1.6, label=method, zorder=2)
            if component == "zz" and row == 0:
                for method, (times, product) in magnetization_products.items():
                    axis.plot(times, product, color=colors[method], linestyle=":",
                              linewidth=1.7, zorder=3)
            if row == 0:
                axis.set_title(rf"$\beta J_Q={beta:g}$")
            if row == 1:
                axis.set_xlabel(r"$t J_Q$")
            axis.set_xlim(0.0, TMAX)
            axis.grid(True, alpha=0.30)
            if column == 0:
                part_label = r"\mathrm{Re}" if part == "real" else r"\mathrm{Im}"
                axis.set_ylabel(rf"${part_label}\,g^{{{component}}}(t)$")

    handles = [Line2D([], [], color=colors[name], linestyle=styles[name], label=name)
               for name in methods]
    if component == "zz":
        handles.append(Line2D([], [], color="black", linestyle=":",
                              label=r"$\langle S^z(t)\rangle\langle S^z(0)\rangle$"))
    axes[0, 0].legend(handles=handles, loc="upper right", fontsize="x-small",
                      framealpha=0.9)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.08)
    case_output = OUTPUT / case_stem
    case_output.mkdir(parents=True, exist_ok=True)
    fig.savefig(case_output / f"g{component}_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    plt.style.use(ROOT / "matplotlibrc")
    benchmark_dir = DATA / "Benchmarks" / "spinDMFT_ana_cont"
    qmc_cont_path = DATA / "Benchmarks" / "QMC_ana_cont" / "realtime_Jminus_hz0.5_N32768_all_betas.hdf5"
    cases = (
        (-2.4495, 0.5, PLOT_BETAS, "JLminus_hz0p5", "realtime_JLminus_hz0.5_all_betas.hdf5",
         qmc_cont_path, None, ("xx", "xy", "zz")),
        (0.0, 0.5, PLOT_BETAS, "JL0_hz0p5", "realtime_JL0_hz0.5_all_betas.hdf5",
         None, 0.5, ("xx", "xy", "zz")),
        (0.0, 0.0, PLOT_BETAS, "JL0_hz0", "realtime_JL0_hz0_all_betas.hdf5",
         None, 0.0, ("xx",)),
    )
    for jl, h_abs, betas, case_stem, spin_file, qmc_file, random_h_abs, components in cases:
        files = staged_files(jl, h_abs, betas)
        spin_cont_path = benchmark_dir / spin_file
        for path in (spin_cont_path, qmc_file):
            if path is not None and not path.is_file():
                raise FileNotFoundError(path)
        for component in components:
            plot_component(component, files, betas, spin_cont_path, qmc_file,
                           random_h_abs, case_stem)


if __name__ == "__main__":
    main()
