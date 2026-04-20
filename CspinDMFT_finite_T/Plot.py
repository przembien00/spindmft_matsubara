from __future__ import annotations

"""Simple generalized plotting script for temperature sweeps.

Edit the user parameters below and run the file. The script scans one folder in
``Data/``, loads matching HDF5 files into nested dictionaries, extracts the
requested correlation, and plots one curve per temperature on the same axis.
"""

import os
from pathlib import Path
import re
from typing import Any

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLBACKEND", "Agg")

import h5py as h5
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


# User parameters
DATA_PROJECT = "SquareLatticeRuns"
SPIN_MODEL = "ISO"  # inserted in the filename as spinmodel=<SPIN_MODEL> before config
CONFIG = "Square2D_2Spin_Bond_NN_AFM_J0p5"
PAIR = "1-2"
CORRELATION_TYPE = "xx"  # xx, xy, xz, yx, yy, yz, zx, zy, zz
DATASET = "auto"  # auto, Re_correlation or Im_correlation
BETA_ARRAY: list[float] | None = None  # set to None to use all available betas
EXTENSION = None  # set to None to accept all matching files
CURVE_COLOR = "limegreen"
PLOT_TAU_OVER_BETA = True
SHOW_ERROR_BARS = False
STACK_ALL_CORRELATIONS = False  # if True, plot all stored components for the symmetry in stacked axes
Y_RANGE = None  # for example (0.1, 0.3)
Y_RANGES = None  # for example {"xx": (0.1, 0.3), "zz": (0.0, 0.2)}
OUTPUT_NAME = None  # set to None for automatic naming
PLOTS_PROJECT = None  # for example "SquareLattice" or "TempSweeps"


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "Data" / DATA_PROJECT
PLOTS_DIR = ROOT / "Plots"
SPINMODEL_RE = re.compile(r"(?:^|__)spinmodel=(.*?)(?=__|$)")
CONFIG_RE = re.compile(r"(?:^|__)config=(.*?)(?=__|$)")
BETA_RE = re.compile(r"(?:^|__)beta=([^_]+)")
MARKERS = ["v", "^", "s", "x", "D", "p", "o", "*", "+", "<", ">"]


def hdf5_to_dict(obj: h5.Group | h5.Dataset) -> Any:
    """Load an HDF5 group or dataset into a nested dictionary."""

    if isinstance(obj, h5.Dataset):
        return obj[...]

    data: dict[str, Any] = {}
    for key, value in obj.items():
        data[key] = hdf5_to_dict(value)

    if obj.attrs:
        data["attrs"] = {key: decode_value(value) for key, value in obj.attrs.items()}

    return data


def decode_value(value: Any) -> Any:
    """Convert HDF5 scalar values to plain Python objects."""

    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        return value.item()
    return value


def load_file(path: Path) -> dict[str, Any]:
    """Load one HDF5 file completely into a dictionary."""

    with h5.File(path, "r") as handle:
        return hdf5_to_dict(handle)


def filename_config(path: Path) -> str | None:
    """Read the configuration label from the filename."""

    match = CONFIG_RE.search(path.stem)
    return match.group(1) if match else None


def filename_spinmodel(path: Path) -> str | None:
    """Read the spin-coupling model from the filename."""

    match = SPINMODEL_RE.search(path.stem)
    return match.group(1) if match else None


def filename_beta(path: Path) -> float | None:
    """Read beta from the filename."""

    match = BETA_RE.search(path.stem)
    return float(match.group(1)) if match else None


def matching_files(
    data_dir: Path,
    spin_model: str | None,
    config: str,
    beta_array: list[float] | None = None,
    extension: str | None = None,
) -> list[Path]:
    """Return all files in the folder that match the selected run setup."""

    files = []
    for path in sorted(data_dir.glob("*.hdf5")):
        if spin_model is not None and filename_spinmodel(path) != spin_model:
            continue
        if filename_config(path) != config:
            continue
        if extension and extension not in path.name:
            continue
        beta = filename_beta(path)
        if beta is None:
            continue
        if beta_array is not None and beta not in beta_array:
            continue
        files.append(path)
    return sorted(files, key=lambda path: filename_beta(path))


def correlation_group(data: dict[str, Any], dataset: str, correlation_type: str) -> dict[str, np.ndarray]:
    """Choose the requested correlation group from the loaded file."""

    return data["results"][resolve_dataset(dataset, correlation_type)]


def std_group(data: dict[str, Any], dataset: str, correlation_type: str) -> dict[str, np.ndarray] | None:
    """Choose the matching standard-deviation group if it exists."""

    runtimedata = data.get("runtimedata", {})
    resolved_dataset = resolve_dataset(dataset, correlation_type)
    if resolved_dataset == "Re_correlation":
        return runtimedata.get("Re_correlation_sample_stds")
    if resolved_dataset == "Im_correlation":
        return runtimedata.get("Im_correlation_sample_stds")
    return None


def resolve_dataset(dataset: str, correlation_type: str) -> str:
    """Resolve the requested dataset, using auto-selection when desired."""

    if dataset != "auto":
        return dataset
    diagonal_components = {"xx", "yy", "zz"}
    if correlation_type in diagonal_components:
        return "Re_correlation"
    return "Im_correlation"


def correlation_labels(symmetry_type: str) -> list[str]:
    """Return the stored component labels for the chosen symmetry scheme."""

    if symmetry_type == "A":
        return ["xx"]
    if symmetry_type == "B":
        return ["xx", "zz"]
    if symmetry_type == "C":
        return ["xx", "xy", "zz"]
    if symmetry_type == "D":
        return ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"]
    raise ValueError(f"unknown correlation symmetry type: {symmetry_type}")


def select_component(values: np.ndarray, correlation_type: str, symmetry_type: str) -> np.ndarray:
    """Extract one correlation component from the stored list."""

    labels = correlation_labels(symmetry_type)
    if values.shape[0] != len(labels):
        raise ValueError(
            f"symmetry type {symmetry_type} expects {len(labels)} components, found {values.shape[0]}"
        )
    if correlation_type not in labels:
        allowed = ", ".join(labels)
        raise ValueError(f"correlation type {correlation_type!r} not available for symmetry {symmetry_type}; use {allowed}")
    return values[labels.index(correlation_type)]


def select_error(std_values: np.ndarray, correlation_type: str, symmetry_type: str, num_samples: float) -> np.ndarray:
    """Convert stored sample standard deviations to standard errors."""

    labels = correlation_labels(symmetry_type)
    if std_values.shape[0] != len(labels):
        raise ValueError(
            f"symmetry type {symmetry_type} expects {len(labels)} components, found {std_values.shape[0]}"
        )
    if correlation_type not in labels:
        allowed = ", ".join(labels)
        raise ValueError(f"correlation type {correlation_type!r} not available for symmetry {symmetry_type}; use {allowed}")
    return std_values[labels.index(correlation_type)] / np.sqrt(num_samples)


def extract_curve(
    data: dict[str, Any],
    pair: str,
    correlation_type: str,
    dataset: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, float]:
    """Extract tau, correlation, optional error bars, and beta from one file."""

    params = data["parameters"]["attrs"]
    beta = float(params["beta"])
    delta_t = float(params["delta_t"])
    num_samples = float(params.get("num_Samples", 1.0))
    symmetry_type = str(params["correlation_symmetry_type"])

    corr = np.asarray(correlation_group(data, dataset, correlation_type)[pair], dtype=float)
    tau = np.arange(corr.shape[1]) * delta_t
    curve = select_component(corr, correlation_type, symmetry_type)

    errors = None
    stds = std_group(data, dataset, correlation_type)
    if stds is not None and pair in stds:
        errors = select_error(np.asarray(stds[pair], dtype=float), correlation_type, symmetry_type, num_samples)

    return tau, curve, errors, beta


def symmetry_type(data: dict[str, Any]) -> str:
    """Read the stored correlation symmetry type from one loaded file."""

    return str(data["parameters"]["attrs"]["correlation_symmetry_type"])


def output_path() -> Path:
    """Choose the output file name."""

    target_dir = PLOTS_DIR
    if PLOTS_PROJECT:
        target_dir = target_dir / PLOTS_PROJECT
    target_dir.mkdir(parents=True, exist_ok=True)

    if OUTPUT_NAME is not None:
        return target_dir / OUTPUT_NAME
    if STACK_ALL_CORRELATIONS:
        filename = f"temperature_sweep__spinmodel={SPIN_MODEL}__config={CONFIG}__pair={PAIR}__stack.png"
    else:
        filename = f"temperature_sweep__spinmodel={SPIN_MODEL}__config={CONFIG}__pair={PAIR}__corr={CORRELATION_TYPE}.png"
    return target_dir / filename


def temperature_label(beta: float) -> str:
    """Build the legend label for one temperature."""

    return rf"$\beta J_Q={beta:.3g}$"

def correlation_ylabel(pair: str, correlation_type: str) -> str:
    """Build the plot y-axis label in correlation notation."""

    label = rf"$g^{{{correlation_type}}}_{{{pair.replace('-', '')}}}(\tau)$"
    if resolve_dataset(DATASET, correlation_type) == "Im_correlation":
        return "Im " + label
    return label


def marker_handles(num_handles: int) -> list[Any]:
    """Create marker-only legend handles."""

    handles = []
    for marker in MARKERS[:num_handles]:
        handle = Line2D([], [], color="black", linestyle="None", marker=marker, markersize=6)
        handles.append(handle)
    return handles


def y_range_for(correlation_type: str) -> tuple[float, float] | None:
    """Return the y-range for one component, if one was configured."""

    if Y_RANGES is not None:
        return Y_RANGES.get(correlation_type)
    return Y_RANGE


def plot_single_axis(files: list[Path]) -> None:
    """Plot one selected correlation component on a single axis."""

    plt.figure(figsize=(8, 3))
    labels = []

    for index, path in enumerate(files):
        data = load_file(path)
        tau, curve, errors, beta = extract_curve(data, PAIR, CORRELATION_TYPE, DATASET)
        x = tau / beta if PLOT_TAU_OVER_BETA else tau
        marker = MARKERS[index % len(MARKERS)]
        labels.append(temperature_label(beta))
        plt.plot(
            x,
            curve,
            color=CURVE_COLOR,
            marker=marker,
            markevery=max(len(x) // 12, 1),
            label=temperature_label(beta),
        )
        if SHOW_ERROR_BARS and errors is not None:
            plt.fill_between(x, curve - errors, curve + errors, alpha=0.15, color=CURVE_COLOR)

    plt.legend(marker_handles(len(files)), labels, fontsize=8)
    plt.xlabel(r"$\tau/\beta$" if PLOT_TAU_OVER_BETA else r"$\tau$")
    plt.ylabel(correlation_ylabel(PAIR, CORRELATION_TYPE))
    plt.grid(alpha=0.25)
    if PLOT_TAU_OVER_BETA:
        plt.xlim(left=0.0, right=1.0)
    axis_range = y_range_for(CORRELATION_TYPE)
    if axis_range is not None:
        plt.ylim(*axis_range)


def plot_stacked_axes(files: list[Path]) -> None:
    """Plot all stored correlation components for one symmetry in stacked axes."""

    first_data = load_file(files[0])
    labels = correlation_labels(symmetry_type(first_data))
    fig, axes = plt.subplots(
        nrows=len(labels),
        ncols=1,
        sharex=True,
        figsize=(8, 3 * len(labels)),
        gridspec_kw={"hspace": 0},
    )

    if len(labels) == 1:
        axes = [axes]

    beta_labels = []
    for index, path in enumerate(files):
        data = load_file(path)
        beta = float(data["parameters"]["attrs"]["beta"])
        marker = MARKERS[index % len(MARKERS)]
        beta_labels.append(temperature_label(beta))

        for axis, corr_label in zip(axes, labels):
            tau, curve, errors, _ = extract_curve(data, PAIR, corr_label, DATASET)
            x = tau / beta if PLOT_TAU_OVER_BETA else tau
            axis.plot(x, curve, color=CURVE_COLOR, marker=marker, markevery=max(len(x) // 12, 1))
            if SHOW_ERROR_BARS and errors is not None:
                axis.fill_between(x, curve - errors, curve + errors, alpha=0.15, color=CURVE_COLOR)

    for axis, corr_label in zip(axes, labels):
        axis.set_ylabel(correlation_ylabel(PAIR, corr_label))
        axis.grid(alpha=0.25)
        if PLOT_TAU_OVER_BETA:
            axis.set_xlim(left=0.0, right=1.0)
        axis_range = y_range_for(corr_label)
        if axis_range is not None:
            axis.set_ylim(*axis_range)

    axes[-1].set_xlabel(r"$\tau/\beta$" if PLOT_TAU_OVER_BETA else r"$\tau$")
    axes[-1].legend(marker_handles(len(files)), beta_labels, fontsize=8, loc=3)

    plt.tight_layout()


def main() -> None:
    """Plot the selected temperature sweep."""

    files = matching_files(DATA_DIR, SPIN_MODEL, CONFIG, BETA_ARRAY, EXTENSION)
    if not files:
        raise FileNotFoundError(f"no matching files found in {DATA_DIR}")

    if STACK_ALL_CORRELATIONS:
        plot_stacked_axes(files)
    else:
        plot_single_axis(files)
    plt.savefig(output_path(), dpi=300)


if __name__ == "__main__":
    main()
