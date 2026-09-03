#!/usr/bin/env python3
"""Compare the physical real-time covariance with an exact doubled embedding."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyze_covariance_frequency import (
    load_data,
    magnitude_db,
    reconstruct_covariance,
)


def physical_real_sector(covariance, num_imaginary, num_real):
    """Return Gamma_RR as [time,(+,spin)/(-,spin),time,(+,spin)/(-,spin)]."""
    result = np.empty((num_real, 6, num_real, 6), dtype=complex)
    branch_offset = (num_imaginary, num_imaginary + num_real)
    for first_branch in range(2):
        first = slice(3 * branch_offset[first_branch],
                      3 * (branch_offset[first_branch] + num_real))
        for second_branch in range(2):
            second = slice(3 * branch_offset[second_branch],
                           3 * (branch_offset[second_branch] + num_real))
            result[:, 3 * first_branch:3 * (first_branch + 1),
                   :, 3 * second_branch:3 * (second_branch + 1)] = (
                covariance[first, second].reshape(num_real, 3, num_real, 3)
            )
    return result


def doubled_circulant_embedding(physical):
    """Embed a block-Toeplitz covariance in a 2N block-circulant covariance.

    Kernel indices 1..N-1 contain positive lags, indices 2N-(1..N-1)
    contain negative lags, and index N is an unused zero separator.  The first
    N time blocks are therefore exactly the requested physical marginal.
    """
    num_real, internal, _, internal_second = physical.shape
    if internal != internal_second:
        raise ValueError("The real-time internal covariance block must be square.")
    embedded_size = 2 * num_real
    kernel = np.zeros((embedded_size, internal, internal), dtype=complex)
    kernel[0] = physical[0, :, 0, :]
    for lag in range(1, num_real):
        kernel[lag] = physical[lag, :, 0, :]
        kernel[embedded_size - lag] = physical[0, :, lag, :]

    row = np.arange(embedded_size)[:, None]
    column = np.arange(embedded_size)[None, :]
    embedded = kernel[(row - column) % embedded_size]
    # Move internal indices between the two time indices.
    embedded = embedded.transpose(0, 2, 1, 3)
    return kernel, embedded


def two_forward_dfts(covariance):
    transformed = np.fft.fft(covariance, axis=0, norm="ortho")
    transformed = np.fft.fft(transformed, axis=2, norm="ortho")
    return np.fft.fftshift(np.fft.fftshift(transformed, axes=0), axes=2)


def paired_frequency_fraction(block, frequencies):
    power = np.sum(np.abs(block) ** 2, axis=(1, 3))
    mask = np.zeros(power.shape, dtype=bool)
    modes = np.fft.fftshift(np.arange(frequencies.size))
    positions = {int(mode): position for position, mode in enumerate(modes)}
    for row, mode in enumerate(modes):
        mask[row, positions[int((-mode) % frequencies.size)]] = True
    return float(power[mask].sum() / power.sum())


def branch_block(transformed, first_branch, second_branch):
    return transformed[
        :, 3 * first_branch:3 * (first_branch + 1),
        :, 3 * second_branch:3 * (second_branch + 1),
    ]


def plot_comparison(output_dir, physical_frequency, embedded_frequency,
                    physical_hat, embedded_hat):
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 6.2), constrained_layout=True)
    image = None
    for row, branch in enumerate((0, 1)):
        for column, (frequencies, transformed, label) in enumerate((
            (physical_frequency, physical_hat, "physical"),
            (embedded_frequency, embedded_hat, "embedded"),
        )):
            block = branch_block(transformed, branch, branch)
            magnitude = np.sqrt(np.sum(np.abs(block) ** 2, axis=(1, 3)))
            image = axes[row, column].pcolormesh(
                frequencies, frequencies,
                magnitude_db(magnitude, np.max(magnitude)),
                shading="auto", cmap="magma", vmin=-100, vmax=0,
            )
            axes[row, column].set_xlim(frequencies[0], frequencies[-1])
            axes[row, column].set_ylim(frequencies[0], frequencies[-1])
            axes[row, column].set_xlabel(r"$\omega'$")
            axes[row, column].set_ylabel(r"$\omega$")
            branch_label = "++" if branch == 0 else "--"
            axes[row, column].text(
                0.03, 0.94, f"{branch_label}, {label}",
                transform=axes[row, column].transAxes, va="top", color="white",
                bbox={"facecolor": "black", "alpha": 0.45,
                      "edgecolor": "none", "pad": 2},
            )
    colorbar = fig.colorbar(image, ax=axes, shrink=0.92)
    colorbar.set_label("magnitude relative to each panel maximum (dB)")
    for extension in ("png", "pdf"):
        fig.savefig(output_dir / f"keldysh_real_circulant_comparison.{extension}", dpi=300)
    plt.close(fig)


def analyze(input_path: Path, output_dir: Path):
    parameters, edge, magnetization, metadata = load_data(input_path)
    covariance, raw_symmetry_error, branch_sizes, _ = reconstruct_covariance(
        parameters, edge, magnetization
    )
    num_imaginary, num_real, _ = branch_sizes
    physical = physical_real_sector(covariance, num_imaginary, num_real)
    kernel, embedded = doubled_circulant_embedding(physical)

    physical_hat = two_forward_dfts(physical)
    embedded_hat = two_forward_dfts(embedded)
    delta_real = float(parameters["delta_real_t"])
    physical_frequency = np.fft.fftshift(
        2 * np.pi * np.fft.fftfreq(num_real, d=delta_real)
    )
    embedded_frequency = np.fft.fftshift(
        2 * np.pi * np.fft.fftfreq(2 * num_real, d=delta_real)
    )

    physical_matrix = physical.transpose(0, 1, 2, 3).reshape(6 * num_real, 6 * num_real)
    embedded_matrix = embedded.reshape(12 * num_real, 12 * num_real)
    physical_marginal = embedded[:num_real, :, :num_real, :]
    marginal_denominator = np.linalg.norm(physical)
    embedded_denominator = np.linalg.norm(embedded_matrix)
    kernel_transpose_error = max(
        np.linalg.norm(kernel[index] - kernel[-index % (2 * num_real)].T)
        for index in range(2 * num_real)
    )

    summary = {
        "input": str(input_path.resolve()),
        "termination": metadata["termination"],
        "physical_real_points": num_real,
        "physical_interval": [0.0, float(parameters["Tmax"])],
        "embedded_real_points": 2 * num_real,
        "embedded_period": 2 * num_real * delta_real,
        "embedded_lag_content": (
            "positive lags 0..T, one zero separator, negative lags -T..-delta_t"
        ),
        "internal_block_dimension": 6,
        "physical_matrix_dimension": int(6 * num_real),
        "embedded_matrix_dimension": int(12 * num_real),
        "reconstructed_raw_symmetry_error": float(raw_symmetry_error),
        "physical_marginal_relative_error": float(
            np.linalg.norm(physical_marginal - physical) / marginal_denominator
        ),
        "embedded_transpose_symmetry_error": float(
            np.linalg.norm(embedded_matrix - embedded_matrix.T) / embedded_denominator
        ),
        "kernel_modular_transpose_max_error": float(kernel_transpose_error),
        "fourier_frobenius_norm_relative_error": float(
            abs(np.sum(np.abs(embedded_hat) ** 2) - np.sum(np.abs(embedded) ** 2))
            / np.sum(np.abs(embedded) ** 2)
        ),
        "paired_frequency_weight": {
            "physical_joint_real_sector": paired_frequency_fraction(
                physical_hat, physical_frequency
            ),
            "embedded_joint_real_sector": paired_frequency_fraction(
                embedded_hat, embedded_frequency
            ),
            "physical_++": paired_frequency_fraction(
                branch_block(physical_hat, 0, 0), physical_frequency
            ),
            "embedded_++": paired_frequency_fraction(
                branch_block(embedded_hat, 0, 0), embedded_frequency
            ),
            "physical_--": paired_frequency_fraction(
                branch_block(physical_hat, 1, 1), physical_frequency
            ),
            "embedded_--": paired_frequency_fraction(
                branch_block(embedded_hat, 1, 1), embedded_frequency
            ),
        },
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    plot_comparison(
        output_dir, physical_frequency, embedded_frequency, physical_hat, embedded_hat
    )
    summary_path = output_dir / "keldysh_real_circulant_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    analyze(args.input, args.output_dir)


if __name__ == "__main__":
    main()
