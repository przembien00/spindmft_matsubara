#!/usr/bin/env python3
"""Fourier-analyze the full Matsubara plus doubled-real Keldysh covariance."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyze_covariance_frequency import load_data, magnitude_db, reconstruct_covariance
from analyze_real_circulant_embedding import (
    doubled_circulant_embedding,
    paired_frequency_fraction,
    physical_real_sector,
    two_forward_dfts,
)


def matsubara_real_block(covariance, num_imaginary, num_real):
    result = np.empty((num_imaginary, 3, num_real, 6), dtype=complex)
    matsubara = slice(0, 3 * num_imaginary)
    for real_branch in range(2):
        start = num_imaginary + real_branch * num_real
        real = slice(3 * start, 3 * (start + num_real))
        result[:, :, :, 3 * real_branch:3 * (real_branch + 1)] = (
            covariance[matsubara, real].reshape(num_imaginary, 3, num_real, 3)
        )
    return result


def kms_extend_mixed_block(physical, embedded_size):
    """Extend t>=0 using X_ab(-t,tau)=X_ba(t,beta-tau)."""
    num_imaginary, _, num_real, internal = physical.shape
    if embedded_size != 2 * num_real or internal != 6:
        raise ValueError("Expected a doubled real grid with two three-spin branches.")
    extended = np.zeros((num_imaginary, 3, embedded_size, 6), dtype=complex)
    extended[:, :, :num_real, :] = physical
    reflection = np.arange(num_imaginary - 1, -1, -1)
    for lag in range(1, num_real):
        for branch in range(2):
            block = physical[reflection, :, lag, 3 * branch:3 * (branch + 1)]
            extended[:, :, embedded_size - lag, 3 * branch:3 * (branch + 1)] = (
                block.transpose(0, 2, 1)
            )
    return extended


def two_grid_dfts(block):
    transformed = np.fft.fft(block, axis=0, norm="ortho")
    transformed = np.fft.fft(transformed, axis=2, norm="ortho")
    return np.fft.fftshift(np.fft.fftshift(transformed, axes=0), axes=2)


def partner_column_order(size):
    modes = np.fft.fftshift(np.arange(size))
    positions = {int(mode): position for position, mode in enumerate(modes)}
    return np.array([positions[int((-mode) % size)] for mode in modes], dtype=int)


def flatten_full_frequency_matrix(mm_hat, mr_hat, rr_hat):
    matsubara_dimension = mm_hat.shape[0] * mm_hat.shape[1]
    real_dimension = rr_hat.shape[0] * rr_hat.shape[1]
    result = np.empty(
        (matsubara_dimension + real_dimension,
         matsubara_dimension + real_dimension), dtype=complex
    )
    mm_matrix = mm_hat.reshape(matsubara_dimension, matsubara_dimension)
    mr_matrix = mr_hat.reshape(matsubara_dimension, real_dimension)
    rr_matrix = rr_hat.reshape(real_dimension, real_dimension)
    result[:matsubara_dimension, :matsubara_dimension] = mm_matrix
    result[:matsubara_dimension, matsubara_dimension:] = mr_matrix
    result[matsubara_dimension:, :matsubara_dimension] = mr_matrix.T
    result[matsubara_dimension:, matsubara_dimension:] = rr_matrix
    return result


def reordered_partner_columns(matrix, num_imaginary, embedded_real):
    num_intervals = num_imaginary - 1
    matsubara_order = partner_column_order(num_intervals)
    real_order = partner_column_order(embedded_real)
    matsubara_columns = (
        matsubara_order[:, None] * 3 + np.arange(3)[None, :]
    ).ravel()
    matsubara_columns = np.concatenate((
        matsubara_columns,
        np.arange(3 * num_intervals, 3 * num_imaginary),
    ))
    real_columns = (
        real_order[:, None] * 6 + np.arange(6)[None, :]
    ).ravel() + 3 * num_imaginary
    return matrix[:, np.concatenate((matsubara_columns, real_columns))]


def endpoint_aware_dfts(matrix, num_imaginary, embedded_real):
    """Transform the N_tau base grid and real grid, leaving beta- untouched."""
    num_intervals = num_imaginary - 1
    matsubara_fft_size = 3 * num_intervals
    matsubara_size = 3 * num_imaginary
    transformed = matrix.copy()
    transformed[:matsubara_fft_size] = np.fft.fft(
        transformed[:matsubara_fft_size].reshape(
            num_intervals, 3, transformed.shape[1]
        ), axis=0, norm="ortho"
    ).reshape(matsubara_fft_size, transformed.shape[1])
    transformed[:, :matsubara_fft_size] = np.fft.fft(
        transformed[:, :matsubara_fft_size].reshape(
            transformed.shape[0], num_intervals, 3
        ), axis=1, norm="ortho"
    ).reshape(transformed.shape[0], matsubara_fft_size)
    transformed[matsubara_size:] = np.fft.fft(
        transformed[matsubara_size:].reshape(
            embedded_real, 6, transformed.shape[1]
        ), axis=0, norm="ortho"
    ).reshape(6 * embedded_real, transformed.shape[1])
    transformed[:, matsubara_size:] = np.fft.fft(
        transformed[:, matsubara_size:].reshape(
            transformed.shape[0], embedded_real, 6
        ), axis=1, norm="ortho"
    ).reshape(transformed.shape[0], 6 * embedded_real)
    matsubara_order = (
        np.fft.fftshift(np.arange(num_intervals))[:, None] * 3
        + np.arange(3)[None, :]
    ).ravel()
    real_order = (
        np.fft.fftshift(np.arange(embedded_real))[:, None] * 6
        + np.arange(6)[None, :]
    ).ravel() + matsubara_size
    order = np.concatenate((
        matsubara_order,
        np.arange(matsubara_fft_size, matsubara_size),
        real_order,
    ))
    return transformed[np.ix_(order, order)]


def plot_full_matrix(output_dir, full_hat, reordered, boundary):
    reference = np.max(np.abs(full_hat))
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), constrained_layout=True)
    image = None
    for ax, matrix, label in zip(
        axes, (full_hat, reordered),
        (r"$F\Gamma F^T$", r"second axis reordered: $q'\mapsto-q'$"),
    ):
        image = ax.imshow(
            magnitude_db(np.abs(matrix), reference), origin="lower",
            interpolation="nearest", aspect="equal", cmap="magma", vmin=-100, vmax=0,
            extent=(-0.5, matrix.shape[1] - 0.5, -0.5, matrix.shape[0] - 0.5),
        )
        ax.axvline(boundary - 0.5, color="white", linewidth=0.6, alpha=0.8)
        ax.axhline(boundary - 0.5, color="white", linewidth=0.6, alpha=0.8)
        ax.set_xlim(-0.5, matrix.shape[1] - 0.5)
        ax.set_ylim(-0.5, matrix.shape[0] - 0.5)
        ax.set_xlabel("frequency-component index")
        ax.set_ylabel("frequency-component index")
        ax.text(
            0.03, 0.96, label, transform=ax.transAxes, va="top", color="white",
            bbox={"facecolor": "black", "alpha": 0.5, "edgecolor": "none", "pad": 2},
        )
        ax.text(
            0.02, 0.02, "M", transform=ax.transAxes, color="white", va="bottom",
        )
        ax.text(
            0.98, 0.02, r"doubled $+/-$", transform=ax.transAxes,
            color="white", ha="right", va="bottom",
        )
    colorbar = fig.colorbar(image, ax=axes, shrink=0.92)
    colorbar.set_label(r"$20\log_{10}(|\widehat\Gamma|/|\widehat\Gamma|_{\max})$ (dB)")
    for extension in ("png", "pdf"):
        fig.savefig(output_dir / f"keldysh_full_matsubara_doubled_fft.{extension}", dpi=300)
    plt.close(fig)


def analyze(input_path: Path, output_dir: Path):
    parameters, edge, magnetization, metadata = load_data(input_path)
    covariance, raw_symmetry_error, branch_sizes, _ = reconstruct_covariance(
        parameters, edge, magnetization, include_beta_endpoint=True
    )
    num_imaginary, num_real, _ = branch_sizes
    embedded_real = 2 * num_real

    mm = covariance[:3 * num_imaginary, :3 * num_imaginary].reshape(
        num_imaginary, 3, num_imaginary, 3
    )
    physical_mr = matsubara_real_block(covariance, num_imaginary, num_real)
    extended_mr = kms_extend_mixed_block(physical_mr, embedded_real)
    physical_rr = physical_real_sector(covariance, num_imaginary, num_real)
    _, extended_rr = doubled_circulant_embedding(physical_rr)

    physical_dimension = 3 * num_imaginary + 6 * num_real
    embedded_dimension = 3 * num_imaginary + 6 * embedded_real
    physical_reordered = np.empty((physical_dimension, physical_dimension), dtype=complex)
    physical_reordered[:3 * num_imaginary, :3 * num_imaginary] = mm.reshape(
        3 * num_imaginary, 3 * num_imaginary
    )
    physical_reordered[:3 * num_imaginary, 3 * num_imaginary:] = physical_mr.reshape(
        3 * num_imaginary, 6 * num_real
    )
    physical_reordered[3 * num_imaginary:, :3 * num_imaginary] = (
        physical_reordered[:3 * num_imaginary, 3 * num_imaginary:].T
    )
    physical_reordered[3 * num_imaginary:, 3 * num_imaginary:] = physical_rr.reshape(
        6 * num_real, 6 * num_real
    )

    enlarged = np.empty((embedded_dimension, embedded_dimension), dtype=complex)
    enlarged[:3 * num_imaginary, :3 * num_imaginary] = mm.reshape(
        3 * num_imaginary, 3 * num_imaginary
    )
    enlarged[:3 * num_imaginary, 3 * num_imaginary:] = extended_mr.reshape(
        3 * num_imaginary, 6 * embedded_real
    )
    enlarged[3 * num_imaginary:, :3 * num_imaginary] = (
        enlarged[:3 * num_imaginary, 3 * num_imaginary:].T
    )
    enlarged[3 * num_imaginary:, 3 * num_imaginary:] = extended_rr.reshape(
        6 * embedded_real, 6 * embedded_real
    )
    physical_indices = np.concatenate((
        np.arange(3 * num_imaginary),
        3 * num_imaginary + np.arange(6 * num_real),
    ))
    physical_marginal = enlarged[np.ix_(physical_indices, physical_indices)]

    full_hat = endpoint_aware_dfts(enlarged, num_imaginary, embedded_real)
    num_intervals = num_imaginary - 1
    mm_base_hat = two_grid_dfts(
        mm[:num_intervals, :, :num_intervals, :]
    )
    rr_hat = two_forward_dfts(extended_rr)
    reordered = reordered_partner_columns(full_hat, num_imaginary, embedded_real)
    mixed_singular_values = np.linalg.svd(
        full_hat[:3 * num_imaginary, 3 * num_imaginary:], compute_uv=False
    )
    mixed_cumulative_weight = (
        np.cumsum(mixed_singular_values**2) / np.sum(mixed_singular_values**2)
    )

    edge_kms_residual = 0.0
    edge_scale = float(np.max(np.abs(edge)))
    direction_partner = (0, 2, 1, 3)
    for direction, partner in enumerate(direction_partner):
        edge_kms_residual = max(
            edge_kms_residual,
            float(np.max(np.abs(edge[0, direction] - edge[0, partner, ::-1]))),
        )
    delta_real = float(parameters["delta_real_t"])
    matsubara_frequency = np.fft.fftshift(
        2 * np.pi * np.fft.fftfreq(
            num_imaginary - 1, d=float(parameters["delta_imag_t"])
        )
    )
    real_frequency = np.fft.fftshift(
        2 * np.pi * np.fft.fftfreq(embedded_real, d=delta_real)
    )

    summary = {
        "input": str(input_path.resolve()),
        "termination": metadata["termination"],
        "mixed_negative_time_extension": "X_ab(-t,tau)=X_ba(t,beta-tau)",
        "edge_grid_KMS_absolute_residual_at_t0": edge_kms_residual,
        "edge_grid_KMS_relative_residual_at_t0": edge_kms_residual / edge_scale,
        "num_matsubara_cells": num_intervals,
        "num_matsubara_one_sided_field_points": num_imaginary,
        "num_physical_real_points": num_real,
        "num_embedded_real_points": embedded_real,
        "embedded_real_period": embedded_real * delta_real,
        "physical_full_dimension": physical_dimension,
        "enlarged_full_dimension": embedded_dimension,
        "reconstructed_raw_symmetry_error": float(raw_symmetry_error),
        "physical_full_marginal_relative_error": float(
            np.linalg.norm(physical_marginal - physical_reordered)
            / np.linalg.norm(physical_reordered)
        ),
        "enlarged_transpose_symmetry_error": float(
            np.linalg.norm(enlarged - enlarged.T) / np.linalg.norm(enlarged)
        ),
        "full_fourier_transpose_symmetry_error": float(
            np.linalg.norm(full_hat - full_hat.T) / np.linalg.norm(full_hat)
        ),
        "fourier_frobenius_norm_relative_error": float(
            abs(np.sum(np.abs(full_hat) ** 2) - np.sum(np.abs(enlarged) ** 2))
            / np.sum(np.abs(enlarged) ** 2)
        ),
        "paired_frequency_weight": {
            "MM_base_grid": paired_frequency_fraction(
                mm_base_hat, matsubara_frequency
            ),
            "doubled_joint_real_sector": paired_frequency_fraction(rr_hat, real_frequency),
        },
        "extended_mixed_singular_weight_rank": {
            str(threshold): int(np.searchsorted(mixed_cumulative_weight, threshold) + 1)
            for threshold in (0.9, 0.99, 0.999, 0.9999, 0.99999)
        },
        "extended_mixed_cumulative_weight_by_rank": {
            str(rank): float(mixed_cumulative_weight[rank - 1])
            for rank in (3, 6, 9, 12) if rank <= mixed_cumulative_weight.size
        },
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    plot_full_matrix(output_dir, full_hat, reordered, 3 * num_imaginary)
    summary_path = output_dir / "keldysh_full_matsubara_doubled_fft_summary.json"
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
