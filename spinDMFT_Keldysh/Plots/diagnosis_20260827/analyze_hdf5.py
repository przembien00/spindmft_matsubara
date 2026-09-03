#!/usr/bin/env python3
"""Read-only diagnostics for the 2026-08-27 spinDMFT_Keldysh data inventory."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import h5py
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "Data"
OUT = Path(__file__).resolve().parent
plt.style.use(ROOT / "matplotlibrc")


def decoded(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode(errors="replace")
    if isinstance(value, np.ndarray) and value.shape == ():
        return decoded(value.item())
    return value


def labels(results: h5py.Group) -> list[str]:
    return [decoded(x) for x in results["correlation_direction_labels"][:]]


def complex_correlations(handle: h5py.File) -> np.ndarray:
    return (
        handle["results/Re_correlation_edge"][:]
        + 1j * handle["results/Im_correlation_edge"][:]
    )


def complex_magnetization(handle: h5py.File) -> np.ndarray:
    return (
        handle["results/magnetization_time_Re"][:]
        + 1j * handle["results/magnetization_time_Im"][:]
    )


def closure_ratio(handle: h5py.File) -> np.ndarray:
    attrs = handle["runtimedata"].attrs
    return np.asarray(attrs["closed_contour_ratio_Re"]) + 1j * np.asarray(
        attrs["closed_contour_ratio_Im"]
    )


def per_time_kms(handle: h5py.File) -> tuple[np.ndarray, np.ndarray]:
    corr = complex_correlations(handle)
    residual = np.conj(corr[:, :, -1]) - corr[:, :, 0]
    magnitude = np.abs(residual)
    component = np.argmax(magnitude, axis=1)
    selected = magnitude[np.arange(magnitude.shape[0]), component]

    re_std = handle["runtimedata/Re_correlation_edge_sample_stds"][:]
    im_std = handle["runtimedata/Im_correlation_edge_sample_stds"][:]
    sigma = np.empty(magnitude.shape[0])
    for t, p in enumerate(component):
        sigma[t] = np.sqrt(
            re_std[t, p, 0] ** 2
            + re_std[t, p, -1] ** 2
            + im_std[t, p, 0] ** 2
            + im_std[t, p, -1] ** 2
        )
    return selected, sigma


def per_time_endpoint_kink(handle: h5py.File) -> np.ndarray:
    corr = complex_correlations(handle)
    # Difference of the last two one-step increments.  This removes the
    # ordinary local slope and isolates the beta-edge curvature/contact.
    kink = corr[:, :, -1] - 2 * corr[:, :, -2] + corr[:, :, -3]
    return np.max(np.abs(kink), axis=1)


def max_with_location(values: np.ndarray) -> tuple[float, tuple[int, ...]]:
    index = np.unravel_index(np.argmax(np.abs(values)), values.shape)
    return float(np.abs(values[index])), index


@dataclass(frozen=True)
class Run:
    label: str
    relative_path: str
    color: str
    linestyle: str = "-"

    @property
    def path(self) -> Path:
        return DATA / self.relative_path


FIELD_COMPARISON = (
    Run(
        r"$h_z=0$",
        "NOT_CONVERGED/spinmodel=ISO__beta=2__factor=fft_wcut=10X.hdf5",
        "blue",
    ),
    Run(
        r"$h_z=0.5$",
        "NOT_CONVERGED/spinmodel=ISO__beta=2__h=z_h_abs=0.5__factor=fft_wcut=10X.hdf5",
        "crimson",
    ),
)

TMAX_COMPARISON = (
    Run(
        r"$T_{\max}=10,\ \Delta t_R=0.05$",
        "spinmodel=ISO__beta=2__h=z_h_abs=0.5__factor=fft_wcut=10__mag=constantX.hdf5",
        "blue",
    ),
    Run(
        r"$T_{\max}=20,\ \Delta t_R=0.10$",
        "NOT_CONVERGED/spinmodel=ISO__beta=2__h=z_h_abs=0.5__factor=fft_wcut=10__mag=constantXX.hdf5",
        "crimson",
    ),
)


def inventory_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in sorted(DATA.rglob("*.hdf5")):
        with h5py.File(path, "r") as handle:
            p = handle["parameters"].attrs
            r = handle["runtimedata"].attrs
            result = handle["results"]
            corr = complex_correlations(handle)
            mag = complex_magnetization(handle)
            closure = closure_ratio(handle)
            kms = np.asarray(r.get("kms_endpoint_residuals", []))
            station = np.asarray(r.get("magnetization_stationarity_residuals", []))
            iteration_error = np.asarray(r.get("absolute_iteration_errors", []))
            approx = np.asarray(r.get("gaussian_covariance_approximation_errors", []))
            phase = np.asarray(r.get("average_phase_magnitudes", []))
            rows.append(
                {
                    "file": str(path.relative_to(DATA)),
                    "schema": (
                        "current_edge_only"
                        if "Re_correlation_center" not in result
                        else "legacy_edge_plus_center"
                    ),
                    "imaginary_field_grid": decoded(p.get("imaginary_field_grid", "")),
                    "equal_time_prescription": decoded(
                        p.get("equal_time_prescription", "not_stored")
                    ),
                    "connected_strategy": decoded(
                        p.get("connected_correlation_strategy", "not_stored")
                    ),
                    "constant_magnetization": decoded(
                        p.get("constant_magnetization_time", "not_stored")
                    ),
                    "labels": ",".join(labels(result)),
                    "beta": p.get("beta", ""),
                    "h_abs": p.get("h_abs", 0.0),
                    "JL": p.get("JL", ""),
                    "JQ": p.get("JQ", ""),
                    "Tmax": p.get("Tmax", ""),
                    "delta_real_t": p.get("delta_real_t", ""),
                    "delta_imag_t": p.get("delta_imag_t", ""),
                    "num_samples": p.get("num_Samples", ""),
                    "iterations": r.get("num_Iterations", ""),
                    "termination": decoded(r.get("termination", "")),
                    "factorization": decoded(p.get("gaussian_factorization", "")),
                    "fft_cutoff": p.get("fft_cross_frequency_cutoff", ""),
                    "iteration_error_last": iteration_error[-1] if iteration_error.size else "",
                    "kms_last": kms[-1] if kms.size else "",
                    "closure_max": np.max(np.abs(closure - 1)),
                    "stationarity_last": station[-1] if station.size else "",
                    "magnetization_drift_max": np.max(np.abs(mag - mag[0])),
                    "phase_last": phase[-1] if phase.size else "",
                    "approximation_error_last": approx[-1] if approx.size else "",
                    "max_Re_correlation_SE": np.max(
                        handle["runtimedata/Re_correlation_edge_sample_stds"][:]
                    ),
                    "edge_shape": "x".join(
                        str(x) for x in result["Re_correlation_edge"].shape
                    ),
                    "finite": bool(
                        np.all(np.isfinite(corr))
                        and np.all(np.isfinite(mag))
                        and np.all(np.isfinite(closure))
                    ),
                }
            )
    return rows


def write_inventory(rows: list[dict[str, object]]) -> None:
    fields = list(rows[0])
    for schema, filename in (
        ("current_edge_only", "current_edge_only_inventory.tsv"),
        ("legacy_edge_plus_center", "legacy_center_grid_inventory.tsv"),
    ):
        with (OUT / filename).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(row for row in rows if row["schema"] == schema)


def selected_metrics(runs: Iterable[Run], filename: str) -> None:
    fields = (
        "label",
        "file",
        "kms_max",
        "kms_t",
        "kms_component",
        "kms_sigma_quadrature",
        "kms_z_quadrature",
        "mag_drift_max",
        "mag_drift_t",
        "mag_component",
        "mag_sigma_quadrature",
        "mag_z_quadrature",
        "closure_max",
        "closure_t",
        "endpoint_kink_max",
        "endpoint_kink_t",
    )
    rows = []
    for run in runs:
        with h5py.File(run.path, "r") as handle:
            times = np.asarray(handle["results"].attrs["real_times"])
            names = labels(handle["results"])
            corr = complex_correlations(handle)
            kms_complex = np.conj(corr[:, :, -1]) - corr[:, :, 0]
            kms_max, (kt, kp) = max_with_location(kms_complex)
            re_std = handle["runtimedata/Re_correlation_edge_sample_stds"][:]
            im_std = handle["runtimedata/Im_correlation_edge_sample_stds"][:]
            kms_sigma = np.sqrt(
                re_std[kt, kp, 0] ** 2
                + re_std[kt, kp, -1] ** 2
                + im_std[kt, kp, 0] ** 2
                + im_std[kt, kp, -1] ** 2
            )
            mag = complex_magnetization(handle)
            delta_mag = mag - mag[0]
            mag_max, (mt, mc) = max_with_location(delta_mag)
            mag_re_std = handle["runtimedata/magnetization_time_Re_sample_stds"][:]
            mag_im_std = handle["runtimedata/magnetization_time_Im_sample_stds"][:]
            mag_sigma = np.sqrt(
                mag_re_std[mt, mc] ** 2
                + mag_re_std[0, mc] ** 2
                + mag_im_std[mt, mc] ** 2
                + mag_im_std[0, mc] ** 2
            )
            closure = closure_ratio(handle)
            ct = int(np.argmax(np.abs(closure - 1)))
            kink = per_time_endpoint_kink(handle)
            et = int(np.argmax(kink))
            rows.append(
                {
                    "label": run.label,
                    "file": run.relative_path,
                    "kms_max": kms_max,
                    "kms_t": times[kt],
                    "kms_component": names[kp],
                    "kms_sigma_quadrature": kms_sigma,
                    "kms_z_quadrature": kms_max / kms_sigma,
                    "mag_drift_max": mag_max,
                    "mag_drift_t": times[mt],
                    "mag_component": "xyz"[mc],
                    "mag_sigma_quadrature": mag_sigma,
                    "mag_z_quadrature": mag_max / mag_sigma,
                    "closure_max": np.abs(closure[ct] - 1),
                    "closure_t": times[ct],
                    "endpoint_kink_max": kink[et],
                    "endpoint_kink_t": times[et],
                }
            )
    with (OUT / filename).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def endpoint_grid_comparison() -> None:
    runs = (
        Run(
            "legacy_cell_centered_dense_16k",
            "spinmodel=ISO__beta=2__h=z_h_abs=0.5X.hdf5",
            "dimgray",
        ),
        Run(
            "current_periodic_edge_fft_cutoff10_80k",
            "NOT_CONVERGED/spinmodel=ISO__beta=2__h=z_h_abs=0.5__factor=fft_wcut=10X.hdf5",
            "crimson",
        ),
    )
    fields = (
        "label",
        "file",
        "component",
        "actual_t",
        "last_increment",
        "endpoint_second_difference",
        "local_second_difference_median",
        "local_second_difference_max",
        "endpoint_over_local_median",
        "endpoint_over_local_max",
    )
    rows = []
    for run in runs:
        with h5py.File(run.path, "r") as handle:
            time = np.asarray(handle["results"].attrs["real_times"])
            t_index = int(np.argmin(np.abs(time - 2.5)))
            corr = complex_correlations(handle)
            for p, component in enumerate(labels(handle["results"])):
                first_difference = np.diff(corr[t_index, p])
                second_difference = np.diff(first_difference)
                local = np.abs(second_difference[-7:-1])
                endpoint = float(np.abs(second_difference[-1]))
                rows.append(
                    {
                        "label": run.label,
                        "file": run.relative_path,
                        "component": component,
                        "actual_t": time[t_index],
                        "last_increment": np.abs(first_difference[-1]),
                        "endpoint_second_difference": endpoint,
                        "local_second_difference_median": np.median(local),
                        "local_second_difference_max": np.max(local),
                        "endpoint_over_local_median": endpoint / np.median(local),
                        "endpoint_over_local_max": endpoint / np.max(local),
                    }
                )
    with (OUT / "endpoint_grid_comparison_t2p5.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def save_figure(
    fig: plt.Figure, stem: str, rect: tuple[float, float, float, float] | None = None
) -> None:
    fig.tight_layout(rect=rect)
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.png", dpi=220)
    plt.close(fig)


def plot_time_diagnostics(runs: Iterable[Run], stem: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.6), sharex=True)
    axes = axes.ravel()
    xmax = 0.0
    for run in runs:
        with h5py.File(run.path, "r") as handle:
            time = np.asarray(handle["results"].attrs["real_times"])
            xmax = max(xmax, float(time[-1]))
            closure = np.abs(closure_ratio(handle) - 1)
            mag = complex_magnetization(handle)
            mag_drift = np.max(np.abs(mag - mag[0]), axis=1)
            kms, kms_sigma = per_time_kms(handle)
            kink = per_time_endpoint_kink(handle)
            axes[0].plot(time, closure, color=run.color, ls=run.linestyle, label=run.label)
            axes[1].plot(time, mag_drift, color=run.color, ls=run.linestyle)
            axes[2].plot(time, kms, color=run.color, ls=run.linestyle)
            axes[2].plot(time, kms_sigma, color=run.color, ls=":", alpha=0.8)
            axes[3].plot(time, kink, color=run.color, ls=run.linestyle)
    axes[0].set_ylabel(r"$|D(t)-1|$")
    axes[1].set_ylabel(r"$\max_a|m_a(t)-m_a(0)|$")
    axes[2].set_ylabel(r"$\max_{ab}|X^{ab}(t,\beta)^*-X^{ab}(t,0)|$")
    axes[3].set_ylabel(r"$\max_{ab}|\Delta_\tau^2X^{ab}|_{\beta}$")
    axes[0].legend(loc="best")
    for axis in axes:
        axis.set_xlim(0, xmax)
        axis.margins(x=0)
    for axis in axes[2:]:
        axis.set_xlabel(r"$tJ$")
    save_figure(fig, stem)


def plot_iteration_history(runs: Iterable[Run], stem: str) -> None:
    attributes = (
        ("absolute_iteration_errors", r"fixed-point error"),
        ("denominator_constancy_residuals", r"$\max_t|D(t)-1|$"),
        ("magnetization_stationarity_residuals", r"$\max_{t,a}|m_a(t)-m_a(0)|$"),
        ("kms_endpoint_residuals", r"KMS endpoint residual"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.5), sharex=True)
    axes = axes.ravel()
    max_iteration = 1
    for run in runs:
        with h5py.File(run.path, "r") as handle:
            attrs = handle["runtimedata"].attrs
            for axis, (name, _) in zip(axes, attributes):
                values = np.asarray(attrs[name])
                iteration = np.arange(1, values.size + 1)
                max_iteration = max(max_iteration, int(iteration[-1]))
                axis.plot(
                    iteration,
                    values,
                    color=run.color,
                    ls=run.linestyle,
                    marker="o",
                    markevery=max(1, values.size // 10),
                    label=run.label,
                )
            threshold = handle["parameters"].attrs["absolute_iteration_error_threshold"]
            axes[0].axhline(threshold, color=run.color, ls=":", alpha=0.8)
    for axis, (_, ylabel) in zip(axes, attributes):
        axis.set_ylabel(ylabel)
        axis.set_xlim(1, max_iteration)
        axis.margins(x=0)
    axes[0].legend(loc="best")
    for axis in axes[2:]:
        axis.set_xlabel("iteration")
    save_figure(fig, stem)


def plot_endpoint_profiles() -> None:
    legacy = Run(
        "cell-centered, dense, 16k",
        "spinmodel=ISO__beta=2__h=z_h_abs=0.5X.hdf5",
        "dimgray",
        "--",
    )
    current = Run(
        "periodic-edge, FFT cutoff 10, 80k",
        "NOT_CONVERGED/spinmodel=ISO__beta=2__h=z_h_abs=0.5__factor=fft_wcut=10X.hdf5",
        "crimson",
        "-",
    )
    components = (("xx", "Re"), ("xy", "Im"), ("yx", "Im"), ("zz", "Re"))
    fig, axes = plt.subplots(2, 2, figsize=(9.1, 6.7), sharex=True)
    axes = axes.ravel()
    for run in (legacy, current):
        with h5py.File(run.path, "r") as handle:
            time = np.asarray(handle["results"].attrs["real_times"])
            t_index = int(np.argmin(np.abs(time - 2.5)))
            beta = float(handle["parameters"].attrs["beta"])
            tau = np.asarray(handle["results"].attrs["imaginary_time_edges"]) / beta
            names = labels(handle["results"])
            re = handle["results/Re_correlation_edge"][:]
            im = handle["results/Im_correlation_edge"][:]
            for axis, (component, part) in zip(axes, components):
                p = names.index(component)
                values = re[t_index, p] if part == "Re" else im[t_index, p]
                axis.plot(tau, values, color=run.color, ls=run.linestyle, label=run.label)
                zoom = getattr(axis, "_diagnostic_inset", None)
                if zoom is None:
                    zoom = inset_axes(axis, width="43%", height="43%", loc="lower left")
                    axis._diagnostic_inset = zoom
                zoom.plot(
                    tau,
                    values,
                    color=run.color,
                    ls=run.linestyle,
                    marker="o",
                    markersize=2.2,
                )
                zoom.set_xlim(0.88, 1.0)
                selected = values[tau >= 0.88]
                padding = 0.08 * max(float(np.ptp(selected)), 1e-12)
                old = getattr(zoom, "_diagnostic_limits", None)
                limits = (float(np.min(selected) - padding), float(np.max(selected) + padding))
                if old is not None:
                    limits = (min(old[0], limits[0]), max(old[1], limits[1]))
                zoom._diagnostic_limits = limits
                zoom.set_ylim(*limits)
                zoom.set_yticklabels([])
                zoom.tick_params(labelsize=7)
    for axis, (component, part) in zip(axes, components):
        axis.set_ylabel(rf"${part}\,X^{{{component}}}(t,\tau)$")
        axis.set_xlim(0, 1)
        axis.margins(x=0)
    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.0),
        ncol=2,
        fontsize=9,
    )
    for axis in axes[2:]:
        axis.set_xlabel(r"$\tau/\beta$")
    save_figure(fig, "endpoint_profiles_legacy_vs_current_t2p5", rect=(0, 0, 1, 0.94))


def main() -> None:
    rows = inventory_rows()
    write_inventory(rows)
    selected_metrics(FIELD_COMPARISON, "matched_current_80k_field_metrics.tsv")
    selected_metrics(TMAX_COMPARISON, "matched_current_80k_Tmax_metrics.tsv")
    endpoint_grid_comparison()
    plot_time_diagnostics(
        FIELD_COMPARISON, "current_80k_field_vs_zero_time_diagnostics"
    )
    plot_iteration_history(
        FIELD_COMPARISON, "current_80k_field_vs_zero_iteration_history"
    )
    plot_time_diagnostics(
        TMAX_COMPARISON, "current_constant_80k_Tmax_time_diagnostics"
    )
    plot_iteration_history(
        TMAX_COMPARISON, "current_constant_80k_Tmax_iteration_history"
    )
    plot_endpoint_profiles()


if __name__ == "__main__":
    main()
