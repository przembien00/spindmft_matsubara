from __future__ import annotations

"""Compare two parameter sets with the same plotting workflow as Plot.py."""

import os
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Any

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLBACKEND", "Agg")

import h5py as h5
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


@dataclass(frozen=True)
class RunSpec:
    data_project: str
    spin_model: str
    config: str
    extension: str | None
    curve_color: str
    curve_style: str
    label: str


# Shared user parameters
PAIR = "1-2"
CORRELATION_TYPE = "xx"  # xx, xy, xz, yx, yy, yz, zx, zy, zz
DATASET = "auto"  # auto, Re_correlation or Im_correlation
BETA_ARRAY: list[float] | None = None  # set to None to use the common available betas
PLOT_TAU_OVER_BETA = True
SHOW_ERROR_BARS = False
STACK_ALL_CORRELATIONS = False
Y_RANGE = None  # for example (0.1, 0.3)
Y_RANGES = None  # for example {"xx": (0.1, 0.3), "zz": (0.0, 0.2)}
OUTPUT_NAME = None
PLOTS_PROJECT = None

# Compared parameter sets
SET_A = RunSpec(
    data_project="SquareLatticeRuns",
    spin_model="ISO",
    config="Square2D_2Spin_Bond_NN_AFM_J0p5",
    extension="cmp",
    curve_color="crimson",
    curve_style="-",
    label="AFM",
)

SET_B = RunSpec(
    data_project="SquareLatticeRuns",
    spin_model="ISO",
    config="Square2D_2Spin_Bond_NN_FM_Jm0p5",
    extension="cmp",
    curve_color="blue",
    curve_style="--",
    label="FM",
)


ROOT = Path(__file__).resolve().parent
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


def matching_files(spec: RunSpec, beta_array: list[float] | None = None) -> list[Path]:
    """Return all files matching one run specification."""

    data_dir = ROOT / "Data" / spec.data_project
    files = []
    for path in sorted(data_dir.glob("*.hdf5")):
        if filename_spinmodel(path) != spec.spin_model:
            continue
        if filename_config(path) != spec.config:
            continue
        if spec.extension and spec.extension not in path.name:
            continue
        beta = filename_beta(path)
        if beta is None:
            continue
        if beta_array is not None and beta not in beta_array:
            continue
        files.append(path)
    return sorted(files, key=lambda path: filename_beta(path))


def file_map(spec: RunSpec, beta_array: list[float] | None = None) -> dict[float, Path]:
    """Map beta values to files for one specification."""

    result: dict[float, Path] = {}
    for path in matching_files(spec, beta_array):
        beta = filename_beta(path)
        if beta is None:
            continue
        result[beta] = path
    return result


def selected_betas(spec_a: RunSpec, spec_b: RunSpec) -> list[float]:
    """Choose the beta values that will be compared."""

    map_a = file_map(spec_a, BETA_ARRAY)
    map_b = file_map(spec_b, BETA_ARRAY)
    common = sorted(set(map_a).intersection(map_b))
    if not common:
        raise FileNotFoundError("no common beta values found for the two compared parameter sets")
    return common


def correlation_group(data: dict[str, Any], dataset: str, correlation_type: str) -> dict[str, np.ndarray]:
    """Choose the requested correlation group from the loaded file."""

    results = data["results"]
    resolved_dataset = resolve_dataset(dataset, correlation_type)
    if resolved_dataset in results:
        return results[resolved_dataset]
    if resolved_dataset == "Re_correlation" and "correlation" in results:
        return results["correlation"]
    raise KeyError(f"results group {resolved_dataset!r} not found")


def std_group(data: dict[str, Any], dataset: str, correlation_type: str) -> dict[str, np.ndarray] | None:
    """Choose the matching standard-deviation group if it exists."""

    runtimedata = data.get("runtimedata", {})
    resolved_dataset = resolve_dataset(dataset, correlation_type)
    if resolved_dataset == "Re_correlation":
        return runtimedata.get("Re_correlation_sample_stds", runtimedata.get("correlation_stds"))
    if resolved_dataset == "Im_correlation":
        return runtimedata.get("Im_correlation_sample_stds")
    return None


def resolve_dataset(dataset: str, correlation_type: str) -> str:
    """Resolve the requested dataset, using auto-selection when desired."""

    if dataset != "auto":
        return dataset
    if correlation_type in {"xx", "yy", "zz"}:
        return "Re_correlation"
    return "Im_correlation"


def correlation_labels(symmetry_type: str) -> list[str]:
    """Return the stored component labels for the chosen symmetry scheme."""

    if symmetry_type == "A":
        return ["xx"]
    if symmetry_type == "B":
        return ["xx", "zz"]
    if symmetry_type == "C":
        return ["xx", "xy", "yx", "zz"]
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
    """Extract the stored standard error for one correlation component."""

    labels = correlation_labels(symmetry_type)
    if std_values.shape[0] != len(labels):
        raise ValueError(
            f"symmetry type {symmetry_type} expects {len(labels)} components, found {std_values.shape[0]}"
        )
    if correlation_type not in labels:
        allowed = ", ".join(labels)
        raise ValueError(f"correlation type {correlation_type!r} not available for symmetry {symmetry_type}; use {allowed}")
    return std_values[labels.index(correlation_type)]


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
    symm = str(params["correlation_symmetry_type"])

    corr = np.asarray(correlation_group(data, dataset, correlation_type)[pair], dtype=float)
    tau = np.arange(corr.shape[1]) * delta_t
    curve = select_component(corr, correlation_type, symm)

    errors = None
    stds = std_group(data, dataset, correlation_type)
    if stds is not None and pair in stds:
        errors = select_error(np.asarray(stds[pair], dtype=float), correlation_type, symm, num_samples)

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
        filename = (
            f"compare__{SET_A.label}__{SET_B.label}__spinmodel={SET_A.spin_model}"
            f"__pair={PAIR}__stack.png"
        )
    else:
        filename = (
            f"compare__{SET_A.label}__{SET_B.label}__spinmodel={SET_A.spin_model}"
            f"__pair={PAIR}__corr={CORRELATION_TYPE}.png"
        )
    return target_dir / filename


def temperature_label(beta: float) -> str:
    """Build the legend label for one temperature."""

    return rf"$\beta J_Q=${beta:.3g}"


def correlation_ylabel(pair: str, correlation_type: str) -> str:
    """Build the plot y-axis label in correlation notation."""

    label = rf"$g^{{{correlation_type}}}_{{{pair.replace('-', '')}}}(\tau)$"
    if resolve_dataset(DATASET, correlation_type) == "Im_correlation":
        return "Im " + label
    return label


def marker_handles(num_handles: int) -> list[Any]:
    """Create marker-only legend handles for beta values."""

    return [Line2D([], [], color="black", linestyle="None", marker=marker, markersize=6) for marker in MARKERS[:num_handles]]


def spec_handles() -> list[Any]:
    """Create line-only legend handles for the two compared parameter sets."""

    return [
        Line2D([], [], color=SET_A.curve_color, linestyle=SET_A.curve_style),
        Line2D([], [], color=SET_B.curve_color, linestyle=SET_B.curve_style),
    ]


def y_range_for(correlation_type: str) -> tuple[float, float] | None:
    """Return the y-range for one component, if one was configured."""

    if Y_RANGES is not None:
        return Y_RANGES.get(correlation_type)
    return Y_RANGE


def plot_one_spec(axis: Any, spec: RunSpec, beta: float, marker: str, correlation_type: str) -> None:
    """Plot one temperature curve for one parameter set on one axis."""

    path = file_map(spec, [beta]).get(beta)
    if path is None:
        raise FileNotFoundError(f"missing beta={beta} for {spec.label}")
    data = load_file(path)
    tau, curve, errors, _ = extract_curve(data, PAIR, correlation_type, DATASET)
    x = tau / beta if PLOT_TAU_OVER_BETA else tau
    axis.plot(
        x,
        curve,
        color=spec.curve_color,
        linestyle=spec.curve_style,
        marker=marker,
        markevery=max(len(x) // 12, 1),
    )
    if SHOW_ERROR_BARS and errors is not None:
        axis.fill_between(x, curve - errors, curve + errors, alpha=0.15, color=spec.curve_color)


def add_legends(axis: Any, betas: list[float]) -> None:
    """Add separate legends for temperatures and compared parameter sets."""

    beta_legend = axis.legend(marker_handles(len(betas)), [temperature_label(beta) for beta in betas], fontsize=8, loc=3)
    axis.add_artist(beta_legend)
    axis.legend(spec_handles(), [SET_A.label, SET_B.label], fontsize=8, loc=4)


def plot_single_axis(betas: list[float]) -> None:
    """Plot one selected correlation component on a single comparison axis."""

    fig, ax = plt.subplots(figsize=(8, 3))
    for index, beta in enumerate(betas):
        marker = MARKERS[index % len(MARKERS)]
        plot_one_spec(ax, SET_A, beta, marker, CORRELATION_TYPE)
        plot_one_spec(ax, SET_B, beta, marker, CORRELATION_TYPE)

    ax.set_xlabel(r"$\tau/\beta$" if PLOT_TAU_OVER_BETA else r"$\tau$")
    ax.set_ylabel(correlation_ylabel(PAIR, CORRELATION_TYPE))
    ax.grid(alpha=0.25)
    if PLOT_TAU_OVER_BETA:
        ax.set_xlim(left=0.0, right=1.0)
    axis_range = y_range_for(CORRELATION_TYPE)
    if axis_range is not None:
        ax.set_ylim(*axis_range)
    add_legends(ax, betas)
    fig.tight_layout()


def plot_stacked_axes(betas: list[float]) -> None:
    """Plot all stored correlation components for one symmetry in stacked axes."""

    first_file = next(iter(file_map(SET_A, [betas[0]]).values()))
    labels = correlation_labels(symmetry_type(load_file(first_file)))
    fig, axes = plt.subplots(
        nrows=len(labels),
        ncols=1,
        sharex=True,
        figsize=(8, 3 * len(labels)),
        gridspec_kw={"hspace": 0},
    )
    if len(labels) == 1:
        axes = [axes]

    for index, beta in enumerate(betas):
        marker = MARKERS[index % len(MARKERS)]
        for axis, corr_label in zip(axes, labels):
            plot_one_spec(axis, SET_A, beta, marker, corr_label)
            plot_one_spec(axis, SET_B, beta, marker, corr_label)

    for axis, corr_label in zip(axes, labels):
        axis.set_ylabel(correlation_ylabel(PAIR, corr_label))
        axis.grid(alpha=0.25)
        if PLOT_TAU_OVER_BETA:
            axis.set_xlim(left=0.0, right=1.0)
        axis_range = y_range_for(corr_label)
        if axis_range is not None:
            axis.set_ylim(*axis_range)

    axes[-1].set_xlabel(r"$\tau/\beta$" if PLOT_TAU_OVER_BETA else r"$\tau$")
    add_legends(axes[-1], betas)
    fig.tight_layout()


def main() -> None:
    """Plot the requested comparison figure."""

    betas = selected_betas(SET_A, SET_B)
    if STACK_ALL_CORRELATIONS:
        plot_stacked_axes(betas)
    else:
        plot_single_axis(betas)
    plt.savefig(output_path(), dpi=300)


if __name__ == "__main__":
    main()
