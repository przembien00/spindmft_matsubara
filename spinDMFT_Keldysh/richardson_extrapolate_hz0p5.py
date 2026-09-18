#!/usr/bin/env python3
"""Richardson-extrapolate the h_z=0.5 real-time magnetization scans."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
BASELINE_DIR = ROOT / "Data" / "Keldysh_stype_C_hz0.5"
DISCRETIZATION_DIR = ROOT / "Data" / "Keldysh_stype_C_hz0.5_discretization"
PLOTS_DIR = ROOT / "Plots"
SAMPLES_PER_CORE = 100_000
EXPECTED_BETAS = (0.2, 0.5, 1.0, 1.5, 2.0, 2.5)
PAIRS = ((200, 400), (300, 600))
BASELINE_FILENAME_RE = re.compile(r"_samples_per_core=100000\.hdf5$")


@dataclass(frozen=True)
class MagnetizationData:
    path: Path
    beta: float
    num_real_steps: int
    num_imag_steps: int
    delta_real_t: float
    delta_imag_t: float
    times: np.ndarray
    values: np.ndarray
    errors: np.ndarray


@dataclass(frozen=True)
class RichardsonData:
    coarse_steps: int
    fine_steps: int
    times: np.ndarray
    values: np.ndarray
    errors: np.ndarray


def _decode(value: object) -> str:
    return value.decode() if isinstance(value, bytes) else str(value)


def read_magnetization(path: Path) -> MagnetizationData:
    with h5py.File(path, "r") as h5file:
        parameters = h5file["parameters"].attrs
        if int(parameters["num_SamplesPerCore"]) != SAMPLES_PER_CORE:
            raise ValueError(f"Unexpected sample count in {path}")
        if int(parameters["num_Cores"]) != 16:
            raise ValueError(f"Unexpected core count in {path}")
        if not np.isclose(float(parameters["Tmax"]), 15.0):
            raise ValueError(f"Unexpected Tmax in {path}")
        if _decode(parameters["sampling_strategy"]) != "independent":
            raise ValueError(f"Unexpected sampling strategy in {path}")
        if _decode(parameters["antithetic_pairs"]) != "no":
            raise ValueError(f"Unexpected antithetic setting in {path}")
        if not _decode(parameters["spin_insertion_strategy"]).startswith("prefix:"):
            raise ValueError(f"Unexpected insertion strategy in {path}")
        if not np.isclose(float(parameters["fft_cross_frequency_cutoff"]), 10.0):
            raise ValueError(f"Unexpected FFT cutoff in {path}")

        beta = float(parameters["beta"])
        num_real_steps = int(parameters["num_RealTimeSteps"])
        num_imag_steps = int(parameters["num_ImagTimeSteps"])
        delta_real_t = float(parameters["delta_real_t"])
        delta_imag_t = float(parameters["delta_imag_t"])
        times = np.asarray(h5file["results"].attrs["real_times"], dtype=float)
        values = np.asarray(h5file["results/Re_magnetization"][:, 2], dtype=float)
        errors = np.asarray(
            h5file["runtimedata/Re_magnetization_sample_stds"][:, 2], dtype=float
        )

    if times.shape != values.shape or values.shape != errors.shape:
        raise ValueError(f"Inconsistent magnetization arrays in {path}")
    return MagnetizationData(
        path,
        beta,
        num_real_steps,
        num_imag_steps,
        delta_real_t,
        delta_imag_t,
        times,
        values,
        errors,
    )


def load_data() -> dict[tuple[float, int], MagnetizationData]:
    paths = list(DISCRETIZATION_DIR.glob("*.hdf5"))
    paths.extend(
        path
        for path in BASELINE_DIR.glob("*.hdf5")
        if BASELINE_FILENAME_RE.search(path.name) is not None
    )
    data: dict[tuple[float, int], MagnetizationData] = {}
    for path in paths:
        item = read_magnetization(path)
        key = (item.beta, item.num_real_steps)
        if key in data:
            raise ValueError(f"Duplicate beta/grid entry: {key}")
        data[key] = item
    return data


def richardson(
    coarse: MagnetizationData, fine: MagnetizationData
) -> RichardsonData:
    if fine.num_real_steps != 2 * coarse.num_real_steps:
        raise ValueError("Richardson pair must halve the real-time step")
    if fine.num_imag_steps != 2 * coarse.num_imag_steps:
        raise ValueError("Richardson pair must halve the imaginary-time step")
    if not np.isclose(fine.delta_real_t, coarse.delta_real_t / 2):
        raise ValueError("Real-time spacings do not form a factor-two pair")
    if not np.isclose(fine.delta_imag_t, coarse.delta_imag_t / 2):
        raise ValueError("Imaginary-time spacings do not form a factor-two pair")
    if not np.allclose(coarse.times, fine.times[::2], rtol=0.0, atol=1e-13):
        raise ValueError("Fine real-time grid does not contain the coarse grid")

    # The runs have independent generated seeds, so no coarse/fine covariance
    # is available for the uncertainty propagation.
    values = (4.0 * fine.values[::2] - coarse.values) / 3.0
    errors = np.sqrt(16.0 * fine.errors[::2] ** 2 + coarse.errors**2) / 3.0
    return RichardsonData(
        coarse.num_real_steps,
        fine.num_real_steps,
        coarse.times,
        values,
        errors,
    )


def common_grid(
    first: RichardsonData, second: RichardsonData
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rounded_first = np.round(first.times, decimals=12)
    rounded_second = np.round(second.times, decimals=12)
    common, first_indices, second_indices = np.intersect1d(
        rounded_first, rounded_second, return_indices=True
    )
    return (
        common,
        first.values[first_indices],
        first.errors[first_indices],
        second.values[second_indices],
        second.errors[second_indices],
    )


def range_and_error(item: RichardsonData) -> tuple[float, float]:
    minimum = int(np.argmin(item.values))
    maximum = int(np.argmax(item.values))
    value = float(item.values[maximum] - item.values[minimum])
    error = float(np.hypot(item.errors[maximum], item.errors[minimum]))
    return value, error


def make_plots(
    extrapolations: dict[float, dict[tuple[int, int], RichardsonData]]
) -> None:
    plt.style.use(ROOT / "matplotlibrc")
    colors = {(200, 400): "C0", (300, 600): "C1"}

    figure, axes = plt.subplots(3, 2, figsize=(7.2, 7.8), sharex=True)
    for axis, beta in zip(axes.flat, EXPECTED_BETAS):
        for pair, item in extrapolations[beta].items():
            label = rf"${pair[0]}\,/\,{pair[1]}$"
            color = colors[pair]
            axis.plot(item.times, item.values, color=color, label=label)
            axis.fill_between(
                item.times,
                item.values - item.errors,
                item.values + item.errors,
                color=color,
                alpha=0.22,
                linewidth=0.0,
            )
        axis.text(0.035, 0.91, rf"$\beta={beta:g}$", transform=axis.transAxes)
        axis.set_xlim(0.0, 15.0)
    axes[0, 0].legend(loc="best")
    for axis in axes[-1, :]:
        axis.set_xlabel(r"$t$")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$\mathrm{Re}\,m^z(t)$")
    figure.tight_layout()
    figure.savefig(
        PLOTS_DIR / "hz0p5_magnetization_z_richardson.png",
        dpi=300,
        bbox_inches="tight",
    )
    figure.savefig(
        PLOTS_DIR / "hz0p5_magnetization_z_richardson.pdf",
        bbox_inches="tight",
    )
    plt.close(figure)

    comparable_betas = [
        beta for beta in EXPECTED_BETAS if len(extrapolations[beta]) == 2
    ]
    figure, axes = plt.subplots(
        len(comparable_betas), 1, figsize=(6.5, 8.0), sharex=True, squeeze=False
    )
    for axis, beta in zip(axes[:, 0], comparable_betas):
        common, first, first_error, second, second_error = common_grid(
            extrapolations[beta][(200, 400)],
            extrapolations[beta][(300, 600)],
        )
        difference = first - second
        difference_error = np.hypot(first_error, second_error)
        axis.axhline(0.0, color="0.25", linewidth=0.8)
        axis.fill_between(
            common,
            -difference_error,
            difference_error,
            color="0.75",
            linewidth=0.0,
        )
        axis.plot(common, difference, color="C3")
        axis.text(0.035, 0.82, rf"$\beta={beta:g}$", transform=axis.transAxes)
        axis.set_xlim(0.0, 15.0)
        axis.set_ylabel(r"$\Delta m^z$")
    axes[-1, 0].set_xlabel(r"$t$")
    figure.tight_layout()
    figure.savefig(
        PLOTS_DIR / "hz0p5_magnetization_z_richardson_pair_difference.png",
        dpi=300,
        bbox_inches="tight",
    )
    figure.savefig(
        PLOTS_DIR / "hz0p5_magnetization_z_richardson_pair_difference.pdf",
        bbox_inches="tight",
    )
    plt.close(figure)


def write_summary(
    extrapolations: dict[float, dict[tuple[int, int], RichardsonData]]
) -> None:
    fields = (
        "beta",
        "pair",
        "range",
        "range_standard_error",
        "pair_difference_rms",
        "pair_difference_max_abs",
        "pair_difference_rms_sigma",
        "pair_difference_max_abs_sigma",
    )
    path = PLOTS_DIR / "hz0p5_magnetization_z_richardson_summary.csv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for beta in EXPECTED_BETAS:
            pair_metrics: dict[str, float | str] = {
                "pair_difference_rms": "",
                "pair_difference_max_abs": "",
                "pair_difference_rms_sigma": "",
                "pair_difference_max_abs_sigma": "",
            }
            if len(extrapolations[beta]) == 2:
                _, first, first_error, second, second_error = common_grid(
                    extrapolations[beta][(200, 400)],
                    extrapolations[beta][(300, 600)],
                )
                difference = first - second
                difference_error = np.hypot(first_error, second_error)
                standardized = difference / difference_error
                pair_metrics = {
                    "pair_difference_rms": float(np.sqrt(np.mean(difference**2))),
                    "pair_difference_max_abs": float(np.max(np.abs(difference))),
                    "pair_difference_rms_sigma": float(
                        np.sqrt(np.mean(standardized**2))
                    ),
                    "pair_difference_max_abs_sigma": float(
                        np.max(np.abs(standardized))
                    ),
                }
            for pair, item in extrapolations[beta].items():
                value, error = range_and_error(item)
                writer.writerow(
                    {
                        "beta": beta,
                        "pair": f"{pair[0]}/{pair[1]}",
                        "range": value,
                        "range_standard_error": error,
                        **pair_metrics,
                    }
                )


def main() -> None:
    PLOTS_DIR.mkdir(exist_ok=True)
    data = load_data()
    extrapolations: dict[float, dict[tuple[int, int], RichardsonData]] = {
        beta: {} for beta in EXPECTED_BETAS
    }
    for beta in EXPECTED_BETAS:
        for pair in PAIRS:
            coarse = data.get((beta, pair[0]))
            fine = data.get((beta, pair[1]))
            if coarse is None or fine is None:
                continue
            extrapolations[beta][pair] = richardson(coarse, fine)
        if not extrapolations[beta]:
            raise ValueError(f"No complete Richardson pair for beta={beta:g}")

    make_plots(extrapolations)
    write_summary(extrapolations)


if __name__ == "__main__":
    main()
