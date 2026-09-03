#!/usr/bin/env python3
"""Construct the Matsubara-bordered real-frequency arrowhead covariance."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyze_covariance_frequency import load_data, magnitude_db, reconstruct_covariance
from analyze_full_matsubara_doubled_fft import (
    flatten_full_frequency_matrix,
    kms_extend_mixed_block,
    matsubara_real_block,
    two_grid_dfts,
)
from analyze_real_circulant_embedding import (
    doubled_circulant_embedding,
    physical_real_sector,
    two_forward_dfts,
)


def frequency_pair_groups(size: int, internal: int, delta: float):
    """Return coordinate groups for {k,-k}, ordered by nonnegative |omega|."""
    shifted_modes = np.fft.fftshift(np.arange(size))
    position = {int(mode): index for index, mode in enumerate(shifted_modes)}
    raw_frequencies = 2 * np.pi * np.fft.fftfreq(size, d=delta)
    groups = []
    visited = set()
    for mode in range(size):
        if mode in visited:
            continue
        partner = (-mode) % size
        visited.update((mode, partner))
        positions = [position[mode]]
        if partner != mode:
            positions.append(position[partner])
        indices = np.concatenate(
            [p * internal + np.arange(internal, dtype=int) for p in positions]
        )
        groups.append(
            {
                "mode": int(mode),
                "partner": int(partner),
                "omega_abs": float(abs(raw_frequencies[mode])),
                "indices": indices,
            }
        )
    groups.sort(key=lambda item: (item["omega_abs"], item["mode"]))
    return groups


def block_diagonal_part(matrix, groups):
    permutation = np.concatenate([group["indices"] for group in groups])
    permuted = matrix[np.ix_(permutation, permutation)]
    result = np.zeros_like(permuted)
    offsets = [0]
    for group in groups:
        offsets.append(offsets[-1] + group["indices"].size)
        block = slice(offsets[-2], offsets[-1])
        result[block, block] = permuted[block, block]
    return permutation, np.asarray(offsets, dtype=int), permuted, result


def schur_complement(A, B_permuted, C_permuted, offsets, rcond):
    schur = A.copy()
    projected_B = np.zeros_like(B_permuted)
    ranks, conditions, inverse_residuals, range_residuals = [], [], [], []
    blocks = []
    for begin, end in zip(offsets[:-1], offsets[1:]):
        C_block = C_permuted[begin:end, begin:end]
        B_block = B_permuted[:, begin:end]
        singular_values = np.linalg.svd(C_block, compute_uv=False)
        cutoff = rcond * singular_values[0] if singular_values.size else 0.0
        rank = int(np.count_nonzero(singular_values > cutoff))
        condition = (
            float(singular_values[0] / singular_values[rank - 1])
            if rank else float("inf")
        )
        inverse = np.linalg.pinv(C_block, rcond=rcond)
        # Preserve the complex-symmetric congruence algebra numerically.
        inverse = 0.5 * (inverse + inverse.T)
        projection = inverse @ C_block
        B_projected = B_block @ projection
        projected_B[:, begin:end] = B_projected
        schur -= B_block @ inverse @ B_block.T
        scale_C = np.linalg.norm(C_block)
        scale_B = np.linalg.norm(B_block)
        inverse_residuals.append(float(
            np.linalg.norm(C_block @ inverse @ C_block - C_block)
            / scale_C if scale_C else 0.0
        ))
        range_residuals.append(float(
            np.linalg.norm(B_projected - B_block) / scale_B if scale_B else 0.0
        ))
        ranks.append(rank)
        conditions.append(condition)
        blocks.append(C_block)
    return {
        "schur": schur,
        "projected_B": projected_B,
        "blocks": blocks,
        "ranks": np.asarray(ranks),
        "conditions": np.asarray(conditions),
        "inverse_residuals": np.asarray(inverse_residuals),
        "range_residuals": np.asarray(range_residuals),
    }


def plot_matrices(output_dir, A, B, C, schur, matsubara_groups, real_groups):
    matsubara_permutation = np.concatenate([g["indices"] for g in matsubara_groups])
    A_plot = A[np.ix_(matsubara_permutation, matsubara_permutation)]
    B_plot = B[matsubara_permutation]
    schur_plot = schur[np.ix_(matsubara_permutation, matsubara_permutation)]
    matrices = (A_plot, B_plot, C, schur_plot)
    labels = (r"$A_M$", r"$B_{MR}$", r"$C_R$", r"$S_M$")
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.3), constrained_layout=True)
    image = None
    for ax, matrix, label in zip(axes.flat, matrices, labels):
        reference = float(np.max(np.abs(matrix)))
        image = ax.imshow(
            magnitude_db(np.abs(matrix), reference),
            origin="lower", interpolation="nearest", aspect="auto",
            cmap="magma", vmin=-100, vmax=0,
            extent=(-0.5, matrix.shape[1] - 0.5, -0.5, matrix.shape[0] - 0.5),
        )
        ax.set_xlim(-0.5, matrix.shape[1] - 0.5)
        ax.set_ylim(-0.5, matrix.shape[0] - 0.5)
        ax.set_xlabel("column index")
        ax.set_ylabel("row index")
        ax.text(
            0.03, 0.95, label, transform=ax.transAxes, va="top", color="white",
            bbox={"facecolor": "black", "alpha": 0.5, "edgecolor": "none", "pad": 2},
        )
    colorbar = fig.colorbar(image, ax=axes, shrink=0.92)
    colorbar.set_label("magnitude relative to each panel maximum (dB)")
    for extension in ("png", "pdf"):
        fig.savefig(output_dir / f"keldysh_arrowhead_matrices.{extension}", dpi=300)
    plt.close(fig)


def plot_mixed_coupling(output_dir, B, matsubara_groups, real_groups):
    values = np.empty((len(matsubara_groups), len(real_groups)))
    for i, matsubara in enumerate(matsubara_groups):
        for j, real in enumerate(real_groups):
            values[i, j] = np.linalg.norm(B[np.ix_(matsubara["indices"], real["indices"])])
    reference = float(np.max(values))
    fig, ax = plt.subplots(figsize=(8.0, 4.3), constrained_layout=True)
    image = ax.imshow(
        magnitude_db(values, reference), origin="lower", interpolation="nearest",
        aspect="auto", cmap="magma", vmin=-100, vmax=0,
        extent=(-0.5, values.shape[1] - 0.5, -0.5, values.shape[0] - 0.5),
    )
    ax.set_xlim(-0.5, values.shape[1] - 0.5)
    ax.set_ylim(-0.5, values.shape[0] - 0.5)
    ax.set_xlabel("real-frequency-pair block")
    ax.set_ylabel("Matsubara-frequency-pair block")
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label(r"$20\log_{10}(\|B_{mn}\|_F/\|B\|_{\max})$ (dB)")
    for extension in ("png", "pdf"):
        fig.savefig(output_dir / f"keldysh_matsubara_real_arrowhead_coupling.{extension}", dpi=300)
    plt.close(fig)
    return values


def analyze(input_path: Path, output_dir: Path, rcond: float):
    parameters, edge, magnetization, metadata = load_data(input_path)
    covariance, raw_symmetry_error, branch_sizes, _ = reconstruct_covariance(
        parameters, edge, magnetization
    )
    num_imaginary, num_real, _ = branch_sizes
    embedded_real = 2 * num_real
    delta_imaginary = float(parameters["delta_imag_t"])
    delta_real = float(parameters["delta_real_t"])

    mm = covariance[:3 * num_imaginary, :3 * num_imaginary].reshape(
        num_imaginary, 3, num_imaginary, 3
    )
    physical_mr = matsubara_real_block(covariance, num_imaginary, num_real)
    extended_mr = kms_extend_mixed_block(physical_mr, embedded_real)
    physical_rr = physical_real_sector(covariance, num_imaginary, num_real)
    _, extended_rr = doubled_circulant_embedding(physical_rr)

    mm_hat = two_grid_dfts(mm)
    mr_hat = two_grid_dfts(extended_mr)
    rr_hat = two_forward_dfts(extended_rr)
    A = mm_hat.reshape(3 * num_imaginary, 3 * num_imaginary)
    B = mr_hat.reshape(3 * num_imaginary, 6 * embedded_real)
    C = rr_hat.reshape(6 * embedded_real, 6 * embedded_real)

    matsubara_groups = frequency_pair_groups(num_imaginary, 3, delta_imaginary)
    real_groups = frequency_pair_groups(embedded_real, 6, delta_real)
    real_permutation, real_offsets, C_permuted, C_block = block_diagonal_part(C, real_groups)
    B_permuted = B[:, real_permutation]
    schur_data = schur_complement(A, B_permuted, C_permuted, real_offsets, rcond)
    schur = schur_data["schur"]

    C_scale = np.linalg.norm(C_permuted)
    offblock_error = float(np.linalg.norm(C_permuted - C_block) / C_scale)
    A_scale = np.linalg.norm(A)
    B_scale = np.linalg.norm(B_permuted)
    schur_scale = np.linalg.norm(schur)
    projected_B_error = float(
        np.linalg.norm(schur_data["projected_B"] - B_permuted) / B_scale
    )
    reconstructed = np.block([
        [A, schur_data["projected_B"]],
        [schur_data["projected_B"].T, C_block],
    ])
    full = flatten_full_frequency_matrix(mm_hat, mr_hat, rr_hat)
    full_real_permutation = np.concatenate((
        np.arange(3 * num_imaginary), 3 * num_imaginary + real_permutation
    ))
    full_permuted = full[np.ix_(full_real_permutation, full_real_permutation)]
    factorization_compatibility_error = float(
        np.linalg.norm(reconstructed - full_permuted) / np.linalg.norm(full_permuted)
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    coupling = plot_mixed_coupling(output_dir, B, matsubara_groups, real_groups)
    plot_matrices(
        output_dir, A, B_permuted, C_block, schur, matsubara_groups, real_groups
    )
    block_sizes = np.diff(real_offsets)
    block_frequencies = np.asarray([g["omega_abs"] for g in real_groups])
    schur_correction = A - schur
    coupling_relative_by_matsubara = coupling / np.max(
        coupling, axis=1, keepdims=True
    )
    np.savez_compressed(
        output_dir / "keldysh_arrowhead_matrices.npz",
        A_M=A,
        B_MR=B_permuted,
        C_R=C_block,
        S_M=schur,
        real_frequency_permutation=real_permutation,
        real_block_offsets=real_offsets,
        real_block_frequencies=block_frequencies,
        real_block_ranks=schur_data["ranks"],
        matsubara_real_block_norms=coupling,
    )

    summary = {
        "input": str(input_path.resolve()),
        "termination": metadata["termination"],
        "pseudoinverse_rcond": rcond,
        "dimensions": {
            "A_M": list(A.shape),
            "B_MR": list(B_permuted.shape),
            "C_R": list(C_block.shape),
            "S_M": list(schur.shape),
        },
        "frequency_pair_blocks": {
            "matsubara_count": len(matsubara_groups),
            "real_count": len(real_groups),
            "real_size_6_count": int(np.count_nonzero(block_sizes == 6)),
            "real_size_12_count": int(np.count_nonzero(block_sizes == 12)),
        },
        "norms": {
            "A_M_frobenius": float(A_scale),
            "B_MR_frobenius": float(B_scale),
            "C_R_frobenius": float(C_scale),
            "S_M_frobenius": float(schur_scale),
            "S_M_over_A_M": float(schur_scale / A_scale),
            "Schur_correction_frobenius": float(np.linalg.norm(schur_correction)),
            "Schur_correction_relative_to_A_M": float(
                np.linalg.norm(schur_correction) / A_scale
            ),
            "Schur_correction_maximum_absolute_entry": float(
                np.max(np.abs(schur_correction))
            ),
        },
        "real_blocks_coupled_per_matsubara_pair_relative_to_its_maximum": {
            str(db): {
                "minimum": int(np.min(counts)),
                "median": float(np.median(counts)),
                "maximum": int(np.max(counts)),
            }
            for db in (-20, -40, -60, -80)
            for counts in [np.sum(
                coupling_relative_by_matsubara >= 10 ** (db / 20), axis=1
            )]
        },
        "real_sector_off_block_relative_norm": offblock_error,
        "real_block_rank": {
            "minimum": int(np.min(schur_data["ranks"])),
            "maximum": int(np.max(schur_data["ranks"])),
            "full_rank_count": int(np.count_nonzero(schur_data["ranks"] == block_sizes)),
        },
        "real_block_condition_number": {
            "minimum": float(np.min(schur_data["conditions"])),
            "median": float(np.median(schur_data["conditions"])),
            "maximum": float(np.max(schur_data["conditions"])),
        },
        "block_pseudoinverse_residual_max": float(
            np.max(schur_data["inverse_residuals"])
        ),
        "mixed_block_range_residual": projected_B_error,
        "arrowhead_schur_factorization_compatibility_error": (
            factorization_compatibility_error
        ),
        "symmetry_errors": {
            "reconstructed_raw_covariance": float(raw_symmetry_error),
            "A_M": float(np.linalg.norm(A - A.T) / A_scale),
            "C_R_block_diagonal": float(np.linalg.norm(C_block - C_block.T) / C_scale),
            "S_M": float(np.linalg.norm(schur - schur.T) / schur_scale),
        },
    }
    (output_dir / "keldysh_arrowhead_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--rcond", type=float, default=1e-12)
    args = parser.parse_args()
    analyze(args.input, args.output_dir, args.rcond)


if __name__ == "__main__":
    main()
