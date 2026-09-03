#!/usr/bin/env python3
"""Reconstruct and Fourier-analyze the spinDMFT_Keldysh pseudo-covariance.

The implementation mirrors Contour/Contour_Kernel.cpp and
Functions/Functions.cpp for the ISO, no-static-noise data used here.  With
V_hat = F V, the unconjugated pseudo-covariance transforms as
Gamma_hat = F Gamma F^T (two forward DFTs, not F Gamma F^dagger).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


BRANCH_NAMES = ("M", "+", "-")
COMPONENTS = "xyz"


def _decoded(value):
    return value.decode() if isinstance(value, bytes) else value


def load_data(path: Path):
    with h5py.File(path, "r") as handle:
        parameters = {key: _decoded(value) for key, value in handle["parameters"].attrs.items()}
        results = handle["results"]
        edge = results["Re_correlation"][...] + 1j * results["Im_correlation"][...]
        labels = [_decoded(value) for value in results["correlation_direction_labels"][...]]
        magnetization = np.asarray(results["Re_magnetization"][0], dtype=float)
        runtime = handle["runtimedata"].attrs
        rank_key = (
            "gaussian_factor_latent_dimensions"
            if "gaussian_factor_latent_dimensions" in runtime
            else "takagi_numerical_ranks"
        )
        metadata = {
            "termination": _decoded(runtime["termination"]),
            "num_iterations": int(runtime["num_Iterations"]),
            "stored_takagi_rank": int(np.asarray(runtime[rank_key])[-1]),
            "stored_covariance_symmetry_error": float(
                np.asarray(runtime["covariance_symmetry_errors"])[-1]
            ),
            "mean_field_distribution_seconds": float(
                handle["timedata"].attrs["mean-field distribution (av)"]
            ),
            "sampling_seconds": float(
                handle["timedata"].attrs["independent complex-field sampling (av)"]
            ),
        }

    if parameters["spin_model"] != "ISO":
        raise ValueError("This reconstruction currently requires spin_model=ISO (D=identity).")
    if parameters["noise"] != "none":
        raise ValueError("This reconstruction currently requires noise=none.")
    if labels != ["xx", "xy", "yx", "zz"]:
        raise ValueError(f"Expected symmetry-C labels [xx,xy,yx,zz], found {labels}.")
    if not (np.isfinite(edge).all() and np.isfinite(magnetization).all()):
        raise ValueError("The selected HDF5 data contain non-finite values.")
    if edge.shape[2] < 2:
        raise ValueError("The imaginary edge grid needs at least one integration interval.")
    return parameters, edge, magnetization, metadata


def reconstruct_covariance(parameters, edge, magnetization,
                           include_beta_endpoint=False):
    num_real = edge.shape[0]
    beta_index = edge.shape[2] - 1
    num_imaginary = edge.shape[2] if include_beta_endpoint else beta_index
    branch_sizes = (num_imaginary, num_real, num_real)
    offsets = np.cumsum((0,) + branch_sizes)
    direction_index = {(0, 0): 0, (1, 1): 0, (0, 1): 1, (1, 0): 2, (2, 2): 3}

    def tensor_value(array, real_time, first_spin, second_spin, imaginary_time):
        direction = direction_index.get((first_spin, second_spin))
        return 0j if direction is None else array[real_time, direction, imaginary_time]

    def edge_value(real_time, first_spin, second_spin, imaginary_edge):
        return tensor_value(edge, real_time, first_spin, second_spin, imaginary_edge)

    def greater_value(first_spin, second_spin, difference):
        if difference >= 0:
            return edge_value(difference, first_spin, second_spin, beta_index)
        return lesser_value(second_spin, first_spin, -difference)

    def lesser_value(first_spin, second_spin, difference):
        if difference >= 0:
            return edge_value(difference, first_spin, second_spin, 0)
        return greater_value(second_spin, first_spin, -difference)

    def branch_correlation(first_branch, first_point, first_spin,
                           second_branch, second_point, second_spin):
        first_matsubara = first_branch == 0
        second_matsubara = second_branch == 0
        if first_matsubara and second_matsubara:
            if first_point > second_point:
                return edge_value(
                    0, second_spin, first_spin, first_point - second_point
                )
            if first_point < second_point:
                return edge_value(
                    0, first_spin, second_spin, second_point - first_point
                )
            return 0.5 * (
                edge_value(0, second_spin, first_spin, 0)
                + edge_value(0, first_spin, second_spin, 0)
            )
        if not first_matsubara and second_matsubara:
            return edge_value(first_point, first_spin, second_spin, second_point)
        if first_matsubara and not second_matsubara:
            return edge_value(second_point, second_spin, first_spin, first_point)

        difference = first_point - second_point
        greater = greater_value(first_spin, second_spin, difference)
        lesser = lesser_value(first_spin, second_spin, difference)
        if first_branch == 1 and second_branch == 1:
            return greater if difference > 0 else lesser if difference < 0 else 0.5 * (greater + lesser)
        if first_branch == 2 and second_branch == 2:
            return lesser if difference > 0 else greater if difference < 0 else 0.5 * (greater + lesser)
        return lesser if (first_branch, second_branch) == (1, 2) else greater

    layout = [
        (branch, point, component)
        for branch, size in enumerate(branch_sizes)
        for point in range(size)
        for component in range(3)
    ]
    jq_squared = float(parameters["JQ"]) ** 2
    raw = np.empty((len(layout), len(layout)), dtype=complex)
    for row, (first_branch, first_point, first_spin) in enumerate(layout):
        for column, (second_branch, second_point, second_spin) in enumerate(layout):
            raw[row, column] = jq_squared * (
                branch_correlation(
                    first_branch, first_point, first_spin,
                    second_branch, second_point, second_spin,
                )
                - magnetization[first_spin] * magnetization[second_spin]
            )

    # Match self_consistent_equations(): construct the upper triangle once and
    # copy it without conjugation into the lower triangle.
    covariance = np.triu(raw) + np.triu(raw, 1).T
    denominator = np.linalg.norm(raw)
    raw_symmetry_error = np.linalg.norm(raw - raw.T) / denominator if denominator else 0.0
    return covariance, raw_symmetry_error, branch_sizes, offsets


def branchwise_fourier(covariance, branch_sizes, offsets, delta_imaginary, delta_real):
    frequencies = tuple(
        np.fft.fftshift(
            2 * np.pi * np.fft.fftfreq(
                size, d=delta_imaginary if branch == 0 else delta_real
            )
        )
        for branch, size in enumerate(branch_sizes)
    )
    blocks = {}
    for first_branch, first_size in enumerate(branch_sizes):
        first_slice = slice(3 * offsets[first_branch], 3 * offsets[first_branch + 1])
        for second_branch, second_size in enumerate(branch_sizes):
            second_slice = slice(3 * offsets[second_branch], 3 * offsets[second_branch + 1])
            block = covariance[first_slice, second_slice].reshape(
                first_size, 3, second_size, 3
            )
            block = np.fft.fft(block, axis=0, norm="ortho")
            block = np.fft.fft(block, axis=2, norm="ortho")
            blocks[first_branch, second_branch] = np.fft.fftshift(
                np.fft.fftshift(block, axes=0), axes=2
            )
    return frequencies, blocks


def paired_frequency_fraction(block, first_frequencies, second_frequencies):
    power = np.sum(np.abs(block) ** 2, axis=(1, 3))
    mask = np.zeros(power.shape, dtype=bool)
    first_modes = np.fft.fftshift(np.arange(first_frequencies.size))
    second_modes = np.fft.fftshift(np.arange(second_frequencies.size))
    second_positions = {int(mode): position for position, mode in enumerate(second_modes)}
    if first_frequencies.size == second_frequencies.size:
        for row, mode in enumerate(first_modes):
            mask[row, second_positions[int((-mode) % second_frequencies.size)]] = True
    else:
        for row, frequency in enumerate(first_frequencies):
            mask[row, np.argmin(np.abs(second_frequencies + frequency))] = True
    return float(power[mask].sum() / power.sum())


def magnitude_db(values, reference, floor=-100.0):
    tiny = np.finfo(float).tiny
    return np.maximum(20 * np.log10(np.maximum(values, tiny) / reference), floor)


def plot_frequency_blocks(output_dir, frequencies, blocks):
    selected = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
    magnitudes = {
        key: np.sqrt(np.sum(np.abs(blocks[key]) ** 2, axis=(1, 3))) for key in selected
    }
    reference = max(np.max(value) for value in magnitudes.values())
    fig, axes = plt.subplots(2, 3, figsize=(10.0, 6.0), constrained_layout=True)
    image = None
    for ax, key in zip(axes.flat, selected):
        first_branch, second_branch = key
        x = frequencies[second_branch]
        y = frequencies[first_branch]
        image = ax.pcolormesh(
            x, y, magnitude_db(magnitudes[key], reference),
            shading="auto", cmap="magma", vmin=-100, vmax=0,
        )
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(y[0], y[-1])
        ax.set_xlabel(r"$\nu'$" if second_branch == 0 else r"$\omega'$" )
        ax.set_ylabel(r"$\nu$" if first_branch == 0 else r"$\omega$" )
        ax.text(
            0.03, 0.94, f"{BRANCH_NAMES[first_branch]}{BRANCH_NAMES[second_branch]}",
            transform=ax.transAxes, va="top", color="white",
            bbox={"facecolor": "black", "alpha": 0.45, "edgecolor": "none", "pad": 2},
        )
    colorbar = fig.colorbar(image, ax=axes, shrink=0.92)
    colorbar.set_label(r"$20\log_{10}(\|\widehat\Gamma_{AB}\|_F/\|\widehat\Gamma\|_{\max})$ (dB)")
    for extension in ("png", "pdf"):
        fig.savefig(output_dir / f"keldysh_covariance_frequency_blocks.{extension}", dpi=300)
    plt.close(fig)


def plot_mixed_components(output_dir, frequencies, blocks):
    mixed = blocks[0, 1]
    selected = ((0, 0), (0, 1), (1, 0), (2, 2))
    reference = max(np.max(np.abs(mixed[:, a, :, b])) for a, b in selected)
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), constrained_layout=True)
    image = None
    for ax, (m_component, real_component) in zip(axes.flat, selected):
        image = ax.pcolormesh(
            frequencies[1], frequencies[0],
            magnitude_db(np.abs(mixed[:, m_component, :, real_component]), reference),
            shading="auto", cmap="magma", vmin=-100, vmax=0,
        )
        ax.set_xlim(frequencies[1][0], frequencies[1][-1])
        ax.set_ylim(frequencies[0][0], frequencies[0][-1])
        ax.set_xlabel(r"$\omega$")
        ax.set_ylabel(r"$\nu$")
        ax.text(
            0.03, 0.94,
            rf"$M_{{{COMPONENTS[m_component]}}}+_{{{COMPONENTS[real_component]}}}$",
            transform=ax.transAxes, va="top", color="white",
            bbox={"facecolor": "black", "alpha": 0.45, "edgecolor": "none", "pad": 2},
        )
    colorbar = fig.colorbar(image, ax=axes, shrink=0.92)
    colorbar.set_label(r"$20\log_{10}(|\widehat\Gamma_{M+}^{ab}|/|\widehat\Gamma_{M+}|_{\max})$ (dB)")
    for extension in ("png", "pdf"):
        fig.savefig(output_dir / f"keldysh_mixed_frequency_components.{extension}", dpi=300)
    plt.close(fig)


def plot_mixed_singular_values(output_dir, singular_values):
    normalized = singular_values / singular_values[0]
    cumulative = np.cumsum(singular_values**2) / np.sum(singular_values**2)
    ranks = np.arange(1, singular_values.size + 1)
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.3), constrained_layout=True)
    axes[0].semilogy(ranks, normalized, color="#0072B2")
    axes[0].set_xlim(ranks[0], ranks[-1])
    axes[0].set_xlabel("singular-value index")
    axes[0].set_ylabel(r"$s_r/s_1$")
    axes[1].plot(ranks, cumulative, color="#D55E00")
    displayed_rank = min(20, ranks[-1])
    axes[1].set_xlim(ranks[0], displayed_rank)
    axes[1].set_ylim(0, 1.005)
    axes[1].set_xlabel("retained rank")
    axes[1].set_ylabel("cumulative squared weight")
    axes[1].axhline(0.99, color="0.4", linestyle="--", linewidth=0.8)
    axes[1].axhline(0.9999, color="0.6", linestyle=":", linewidth=0.8)
    for threshold, color in ((0.99, "#0072B2"), (0.9999, "#009E73")):
        retained_rank = np.searchsorted(cumulative, threshold) + 1
        axes[1].plot(retained_rank, cumulative[retained_rank - 1], "o", color=color)
    for extension in ("png", "pdf"):
        fig.savefig(output_dir / f"keldysh_mixed_singular_values.{extension}", dpi=300)
    plt.close(fig)


def analyze(input_path: Path, output_dir: Path):
    parameters, edge, magnetization, metadata = load_data(input_path)
    covariance, raw_symmetry_error, branch_sizes, offsets = reconstruct_covariance(
        parameters, edge, magnetization
    )
    frequencies, blocks = branchwise_fourier(
        covariance, branch_sizes, offsets,
        float(parameters["delta_imag_t"]), float(parameters["delta_real_t"]),
    )

    mixed = blocks[0, 1].reshape(3 * branch_sizes[0], 3 * branch_sizes[1])
    singular_values = np.linalg.svd(mixed, compute_uv=False)
    cumulative = np.cumsum(singular_values**2) / np.sum(singular_values**2)
    mixed_power = np.sum(np.abs(blocks[0, 1]) ** 2, axis=(1, 3))
    sorted_power = np.sort(mixed_power.ravel())[::-1]
    cumulative_bins = np.cumsum(sorted_power) / np.sum(sorted_power)
    maximum_index = np.unravel_index(np.argmax(mixed_power), mixed_power.shape)
    transformed_norm_squared = sum(
        float(np.sum(np.abs(block) ** 2)) for block in blocks.values()
    )
    covariance_norm_squared = float(np.sum(np.abs(covariance) ** 2))

    summary = {
        "input": str(input_path.resolve()),
        "termination": metadata["termination"],
        "num_iterations": metadata["num_iterations"],
        "num_matsubara_edge_points": branch_sizes[0],
        "num_real_points": branch_sizes[1],
        "covariance_dimension": int(covariance.shape[0]),
        "stored_takagi_rank": metadata["stored_takagi_rank"],
        "stored_covariance_symmetry_error": metadata["stored_covariance_symmetry_error"],
        "reconstructed_raw_symmetry_error": float(raw_symmetry_error),
        "mean_field_distribution_seconds": metadata["mean_field_distribution_seconds"],
        "sampling_seconds": metadata["sampling_seconds"],
        "fourier_convention": "Gamma_hat = F Gamma F^T with orthonormal forward DFTs",
        "fourier_frobenius_norm_relative_error": float(
            abs(transformed_norm_squared - covariance_norm_squared) / covariance_norm_squared
        ),
        "paired_frequency_weight": {
            "MM": paired_frequency_fraction(blocks[0, 0], frequencies[0], frequencies[0]),
            "++": paired_frequency_fraction(blocks[1, 1], frequencies[1], frequencies[1]),
            "+-": paired_frequency_fraction(blocks[1, 2], frequencies[1], frequencies[2]),
            "--": paired_frequency_fraction(blocks[2, 2], frequencies[2], frequencies[2]),
        },
        "mixed_Mplus_Mminus_max_difference": float(np.max(np.abs(blocks[0, 1] - blocks[0, 2]))),
        "mixed_singular_weight_rank": {
            str(threshold): int(np.searchsorted(cumulative, threshold) + 1)
            for threshold in (0.9, 0.95, 0.99, 0.999, 0.9999)
        },
        "mixed_cumulative_weight_by_rank": {
            str(rank): float(cumulative[rank - 1]) for rank in (1, 2, 3, 6)
        },
        "mixed_relative_singular_value_rank": {
            str(tolerance): int(np.count_nonzero(singular_values / singular_values[0] > tolerance))
            for tolerance in (1e-2, 1e-3, 1e-4, 1e-6, 1e-8)
        },
        "mixed_top_bin_weight": {
            str(fraction): float(cumulative_bins[max(1, int(np.ceil(fraction * sorted_power.size))) - 1])
            for fraction in (0.001, 0.005, 0.01, 0.05, 0.1)
        },
        "mixed_maximum_frequency": {
            "matsubara_angular_frequency": float(frequencies[0][maximum_index[0]]),
            "real_angular_frequency": float(frequencies[1][maximum_index[1]]),
        },
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    plot_frequency_blocks(output_dir, frequencies, blocks)
    plot_mixed_components(output_dir, frequencies, blocks)
    plot_mixed_singular_values(output_dir, singular_values)
    summary_path = output_dir / "keldysh_covariance_frequency_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="spinDMFT_Keldysh HDF5 result")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    analyze(args.input, args.output_dir)


if __name__ == "__main__":
    main()
