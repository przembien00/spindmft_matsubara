#!/usr/bin/env python3
"""Construct the eta/nu triangular factor and compare it with dense Takagi."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyze_arrowhead_schur import block_diagonal_part, frequency_pair_groups
from analyze_covariance_frequency import load_data, magnitude_db, reconstruct_covariance
from analyze_full_matsubara_doubled_fft import (
    kms_extend_mixed_block,
    matsubara_real_block,
    two_grid_dfts,
)
from analyze_real_circulant_embedding import (
    doubled_circulant_embedding,
    physical_real_sector,
    two_forward_dfts,
)


def takagi_factor(covariance):
    """Mirror the solver's positive-eigenvalue real-lift Takagi construction."""
    size = covariance.shape[0]
    lift = np.block([
        [covariance.real, covariance.imag],
        [covariance.imag, -covariance.real],
    ])
    eigenvalues, eigenvectors = np.linalg.eigh(lift)
    sigma_max = max(0.0, float(eigenvalues[-1])) if eigenvalues.size else 0.0
    tolerance = (
        np.finfo(float).eps * 100.0 * 2 * size * max(1.0, sigma_max)
    )
    selected = eigenvalues > tolerance
    vectors = eigenvectors[:, selected]
    factor = (
        vectors[:size] + 1j * vectors[size:]
    ) * np.sqrt(eigenvalues[selected])[None, :]
    return factor, tolerance


def eta_nu_transform(number_of_frequencies):
    identity = np.eye(3)
    single = np.block([[0.5 * identity, 0.5 * identity],
                       [identity, -identity]])
    interleaved = np.kron(np.eye(number_of_frequencies), single)
    eta = np.concatenate([
        np.arange(6 * f, 6 * f + 3) for f in range(number_of_frequencies)
    ])
    nu = np.concatenate([
        np.arange(6 * f + 3, 6 * f + 6) for f in range(number_of_frequencies)
    ])
    order = np.concatenate((eta, nu))
    return interleaved[order]


def build_frequency_factor(A, B, C, real_groups, real_permutation, offsets):
    matsubara_factor, matsubara_tolerance = takagi_factor(A)
    block_data = []
    total_columns = matsubara_factor.shape[1]
    for group, begin, end in zip(real_groups, offsets[:-1], offsets[1:]):
        old_C = C[begin:end, begin:end]
        old_B = B[:, begin:end]
        number_of_frequencies = group["indices"].size // 6
        transform = eta_nu_transform(number_of_frequencies)
        transformed_C = transform @ old_C @ transform.T
        transformed_B = old_B @ transform.T
        half = transformed_C.shape[0] // 2
        S = transformed_C[:half, :half]
        R = transformed_C[:half, half:]
        M_nu = transformed_B[:, half:]
        nu_nu = transformed_C[half:, half:]
        U, singular_values, Vh = np.linalg.svd(R, full_matrices=False)
        tolerance = (
            np.finfo(float).eps * 100.0 * max(R.shape)
            * max(1.0, float(singular_values[0]))
        )
        rank = int(np.count_nonzero(singular_values > tolerance))
        root = np.sqrt(singular_values[:rank])
        P = U[:, :rank] * root[None, :]
        Q = Vh[:rank].T * root[None, :]
        block_data.append({
            "begin": int(begin), "end": int(end), "transform": transform,
            "S": S, "R": R, "B": transformed_B[:, :half],
            "M_nu_relative": float(
                np.linalg.norm(M_nu) / np.linalg.norm(transformed_B)
                if np.linalg.norm(transformed_B) else 0.0
            ),
            "nu_nu_relative": float(
                np.linalg.norm(nu_nu) / np.linalg.norm(transformed_C)
                if np.linalg.norm(transformed_C) else 0.0
            ),
            "P": P, "Q": Q, "rank": rank,
            "condition": float(singular_values[0] / singular_values[rank - 1]),
            "column_begin": total_columns,
            "column_end": total_columns + 2 * rank,
        })
        total_columns += 2 * rank

    matsubara_size = A.shape[0]
    real_size = C.shape[0]
    factor_frequency_permuted = np.zeros(
        (matsubara_size + real_size, total_columns), dtype=complex
    )
    factor_frequency_permuted[:matsubara_size, :matsubara_factor.shape[1]] = (
        matsubara_factor
    )
    root_two = np.sqrt(2.0)
    for item in block_data:
        begin, end = item["begin"], item["end"]
        columns = slice(item["column_begin"], item["column_end"])
        rank = item["rank"]
        eta_zero = np.hstack((item["P"] / root_two,
                              1j * item["P"] / root_two))
        nu = np.hstack((item["Q"] / root_two,
                        -1j * item["Q"] / root_two))
        response_solve = np.linalg.solve(item["R"].T, nu)
        eta = eta_zero + 0.5 * item["S"] @ response_solve
        factor_frequency_permuted[:matsubara_size, columns] += (
            item["B"] @ response_solve
        )
        transformed_rows = np.vstack((eta, nu))
        old_rows = np.linalg.solve(item["transform"], transformed_rows)
        factor_frequency_permuted[
            matsubara_size + begin:matsubara_size + end, columns
        ] = old_rows

    factor_frequency = np.empty_like(factor_frequency_permuted)
    factor_frequency[:matsubara_size] = factor_frequency_permuted[:matsubara_size]
    real_inverse_permutation = np.empty_like(real_permutation)
    real_inverse_permutation[real_permutation] = np.arange(real_permutation.size)
    factor_frequency[matsubara_size:] = factor_frequency_permuted[
        matsubara_size + real_inverse_permutation
    ]
    return factor_frequency, block_data, matsubara_tolerance


def frequency_to_physical_factor(factor, num_imaginary, num_real, embedded_real):
    matsubara_size = 3 * num_imaginary
    matsubara_frequency = factor[:matsubara_size].reshape(
        num_imaginary, 3, factor.shape[1]
    )
    matsubara_time = np.fft.ifft(
        np.fft.ifftshift(matsubara_frequency, axes=0), axis=0, norm="ortho"
    )
    real_frequency = factor[matsubara_size:].reshape(
        embedded_real, 6, factor.shape[1]
    )
    real_time = np.fft.ifft(
        np.fft.ifftshift(real_frequency, axes=0), axis=0, norm="ortho"
    )[:num_real]
    return np.vstack((
        matsubara_time.reshape(3 * num_imaginary, factor.shape[1]),
        real_time[:, :3].reshape(3 * num_real, factor.shape[1]),
        real_time[:, 3:].reshape(3 * num_real, factor.shape[1]),
    ))


def dense_hermitian_covariance(covariance):
    square = covariance @ covariance.conj().T
    eigenvalues, eigenvectors = np.linalg.eigh(square)
    singular_values = np.sqrt(np.maximum(eigenvalues, 0.0))
    hermitian = (eigenvectors * singular_values[None, :]) @ eigenvectors.conj().T
    return hermitian, singular_values


def balanced_cross_columns(cross_covariance):
    U, singular_values, Vh = np.linalg.svd(cross_covariance, full_matrices=False)
    tolerance = (
        np.finfo(float).eps * 100.0 * max(cross_covariance.shape)
        * max(1.0, float(singular_values[0]))
    )
    rank = int(np.count_nonzero(singular_values > tolerance))
    root = np.sqrt(singular_values[:rank])
    P = U[:, :rank] * root[None, :]
    Q = Vh[:rank].T * root[None, :]
    return np.hstack((P / np.sqrt(2.0), 1j * P / np.sqrt(2.0))), np.hstack((
        Q / np.sqrt(2.0), -1j * Q / np.sqrt(2.0)
    ))


def additive_pq_factor(covariance, num_imaginary, num_real):
    """Reproduce the solver's separate A/S Takagi plus balanced PQ factors."""
    matsubara = 3 * num_imaginary
    real = 3 * num_real
    plus = matsubara
    minus = matsubara + real
    A = covariance[:matsubara, :matsubara]
    M_plus = covariance[:matsubara, plus:plus + real]
    M_minus = covariance[:matsubara, minus:minus + real]
    pp = covariance[plus:plus + real, plus:plus + real]
    pm = covariance[plus:plus + real, minus:minus + real]
    mp = covariance[minus:minus + real, plus:plus + real]
    mm = covariance[minus:minus + real, minus:minus + real]
    B = 0.5 * (M_plus + M_minus)
    S = 0.25 * (pp + pm + mp + mm)
    R = 0.5 * (pp - pm + mp - mm)
    L_A, _ = takagi_factor(A)
    L_S, _ = takagi_factor(S)
    M_B, eta_B = balanced_cross_columns(B)
    eta_R, nu_R = balanced_cross_columns(R)
    columns = L_A.shape[1] + L_S.shape[1] + M_B.shape[1] + eta_R.shape[1]
    transformed = np.zeros((matsubara + 2 * real, columns), dtype=complex)
    cursor = 0
    transformed[:matsubara, cursor:cursor + L_A.shape[1]] = L_A
    cursor += L_A.shape[1]
    transformed[matsubara:matsubara + real, cursor:cursor + L_S.shape[1]] = L_S
    cursor += L_S.shape[1]
    transformed[:matsubara, cursor:cursor + M_B.shape[1]] = M_B
    transformed[matsubara:matsubara + real, cursor:cursor + M_B.shape[1]] = eta_B
    cursor += M_B.shape[1]
    transformed[matsubara:matsubara + real, cursor:cursor + eta_R.shape[1]] = eta_R
    transformed[matsubara + real:, cursor:cursor + eta_R.shape[1]] = nu_R
    result = np.empty_like(transformed)
    result[:matsubara] = transformed[:matsubara]
    eta = transformed[matsubara:matsubara + real]
    nu = transformed[matsubara + real:]
    result[plus:plus + real] = eta + 0.5 * nu
    result[minus:minus + real] = eta - 0.5 * nu
    return result


def retained_rank(eigenvalues, fraction):
    ordered = np.sort(np.maximum(eigenvalues, 0.0))[::-1]
    cumulative = np.cumsum(ordered)
    return int(np.searchsorted(cumulative, fraction * cumulative[-1]) + 1)


def plot_comparison(output_dir, dense_H, structured_H, pq_H, dense_eigenvalues,
                    structured_eigenvalues, pq_eigenvalues, boundaries):
    dense_ordered = np.sort(dense_eigenvalues)[::-1]
    structured_ordered = np.sort(structured_eigenvalues)[::-1]
    pq_ordered = np.sort(pq_eigenvalues)[::-1]
    index = np.arange(1, dense_ordered.size + 1)
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.0), constrained_layout=True)

    axes[0, 0].semilogy(index, dense_ordered, color="#0072B2", label="dense Takagi")
    axes[0, 0].semilogy(index, structured_ordered, color="#D55E00",
                        linestyle="--", label=r"triangular $\eta,\nu$")
    axes[0, 0].semilogy(index, pq_ordered, color="#009E73",
                        linestyle=":", label="additive PQ")
    axes[0, 0].set_xlim(index[0], index[-1])
    axes[0, 0].set_xlabel("Hermitian-variance eigenvalue index")
    axes[0, 0].set_ylabel(r"eigenvalue of $LL^\dagger$")
    axes[0, 0].legend(frameon=False)

    dense_cumulative = np.cumsum(dense_ordered) / np.sum(dense_ordered)
    structured_cumulative = np.cumsum(structured_ordered) / np.sum(structured_ordered)
    pq_cumulative = np.cumsum(pq_ordered) / np.sum(pq_ordered)
    axes[0, 1].plot(index, dense_cumulative, color="#0072B2")
    axes[0, 1].plot(index, structured_cumulative, color="#D55E00", linestyle="--")
    axes[0, 1].plot(index, pq_cumulative, color="#009E73", linestyle=":")
    axes[0, 1].set_xlim(index[0], index[-1])
    axes[0, 1].set_ylim(0, 1.005)
    axes[0, 1].set_xlabel("retained eigenmodes")
    axes[0, 1].set_ylabel("cumulative Hermitian variance")

    diagonal_dense = np.real(np.diag(dense_H))
    diagonal_structured = np.real(np.diag(structured_H))
    diagonal_pq = np.real(np.diag(pq_H))
    field_index = np.arange(diagonal_dense.size)
    axes[1, 0].semilogy(field_index, diagonal_dense, color="#0072B2")
    axes[1, 0].semilogy(field_index, diagonal_structured, color="#D55E00",
                        linestyle="--")
    axes[1, 0].semilogy(field_index, diagonal_pq, color="#009E73", linestyle=":")
    for boundary in boundaries:
        axes[1, 0].axvline(boundary - 0.5, color="0.55", linewidth=0.7)
    axes[1, 0].set_xlim(field_index[0], field_index[-1])
    axes[1, 0].set_xlabel(r"physical field index $(M,+,-)$")
    axes[1, 0].set_ylabel(r"diagonal of $LL^\dagger$")

    difference = np.abs(structured_H - dense_H)
    reference = float(np.max(np.abs(dense_H)))
    image = axes[1, 1].imshow(
        magnitude_db(difference, reference), origin="lower", interpolation="nearest",
        aspect="equal", cmap="magma", vmin=-100, vmax=0,
        extent=(-0.5, difference.shape[1] - 0.5,
                -0.5, difference.shape[0] - 0.5),
    )
    for boundary in boundaries:
        axes[1, 1].axvline(boundary - 0.5, color="white", linewidth=0.5)
        axes[1, 1].axhline(boundary - 0.5, color="white", linewidth=0.5)
    axes[1, 1].set_xlim(-0.5, difference.shape[1] - 0.5)
    axes[1, 1].set_ylim(-0.5, difference.shape[0] - 0.5)
    axes[1, 1].set_xlabel("physical field index")
    axes[1, 1].set_ylabel("physical field index")
    axes[1, 1].text(
        0.03, 0.95, r"$|H_{\eta\nu}-H_{\rm dense}|$",
        transform=axes[1, 1].transAxes, va="top", color="white",
        bbox={"facecolor": "black", "alpha": 0.5, "edgecolor": "none", "pad": 2},
    )
    colorbar = fig.colorbar(image, ax=axes[1, 1], shrink=0.9)
    colorbar.set_label("difference relative to dense maximum (dB)")
    for extension in ("png", "pdf"):
        fig.savefig(output_dir / f"keldysh_eta_nu_factor_comparison.{extension}", dpi=300)
    plt.close(fig)


def analyze(input_path: Path, output_dir: Path):
    parameters, edge, magnetization, metadata = load_data(input_path)
    covariance, _, branch_sizes, _ = reconstruct_covariance(
        parameters, edge, magnetization
    )
    num_imaginary, num_real, _ = branch_sizes
    embedded_real = 2 * num_real
    delta_real = float(parameters["delta_real_t"])

    mm = covariance[:3 * num_imaginary, :3 * num_imaginary].reshape(
        num_imaginary, 3, num_imaginary, 3
    )
    physical_mr = matsubara_real_block(covariance, num_imaginary, num_real)
    extended_mr = kms_extend_mixed_block(physical_mr, embedded_real)
    physical_rr = physical_real_sector(covariance, num_imaginary, num_real)
    _, extended_rr = doubled_circulant_embedding(physical_rr)
    A = two_grid_dfts(mm).reshape(3 * num_imaginary, 3 * num_imaginary)
    B = two_grid_dfts(extended_mr).reshape(3 * num_imaginary, 6 * embedded_real)
    C = two_forward_dfts(extended_rr).reshape(6 * embedded_real, 6 * embedded_real)

    real_groups = frequency_pair_groups(embedded_real, 6, delta_real)
    real_permutation, offsets, C_permuted, _ = block_diagonal_part(C, real_groups)
    B_permuted = B[:, real_permutation]
    frequency_factor, block_data, matsubara_tolerance = build_frequency_factor(
        A, B_permuted, C_permuted, real_groups, real_permutation, offsets
    )
    physical_factor = frequency_to_physical_factor(
        frequency_factor, num_imaginary, num_real, embedded_real
    )
    reconstructed = physical_factor @ physical_factor.T
    pseudo_error = float(
        np.linalg.norm(reconstructed - covariance) / np.linalg.norm(covariance)
    )

    structured_H = physical_factor @ physical_factor.conj().T
    pq_factor = additive_pq_factor(covariance, num_imaginary, num_real)
    pq_pseudo_error = float(
        np.linalg.norm(pq_factor @ pq_factor.T - covariance) / np.linalg.norm(covariance)
    )
    pq_H = pq_factor @ pq_factor.conj().T
    dense_H, dense_eigenvalues = dense_hermitian_covariance(covariance)
    structured_eigenvalues = np.maximum(np.linalg.eigvalsh(structured_H), 0.0)
    pq_eigenvalues = np.maximum(np.linalg.eigvalsh(pq_H), 0.0)
    dense_trace = float(np.trace(dense_H).real)
    structured_trace = float(np.trace(structured_H).real)
    dense_max = float(np.max(dense_eigenvalues))
    structured_max = float(np.max(structured_eigenvalues))
    pq_trace = float(np.trace(pq_H).real)
    pq_max = float(np.max(pq_eigenvalues))
    difference = np.linalg.norm(structured_H - dense_H) / np.linalg.norm(dense_H)
    diagonal_ratio = np.real(np.diag(structured_H)) / np.real(np.diag(dense_H))
    pq_difference = np.linalg.norm(pq_H - dense_H) / np.linalg.norm(dense_H)
    pq_diagonal_ratio = np.real(np.diag(pq_H)) / np.real(np.diag(dense_H))

    output_dir.mkdir(parents=True, exist_ok=True)
    boundaries = (3 * num_imaginary, 3 * num_imaginary + 3 * num_real)
    plot_comparison(
        output_dir, dense_H, structured_H, pq_H, dense_eigenvalues,
        structured_eigenvalues, pq_eigenvalues, boundaries
    )
    np.savez_compressed(
        output_dir / "keldysh_eta_nu_triangular_factor.npz",
        L_physical=physical_factor,
        H_dense=dense_H,
        H_eta_nu=structured_H,
        H_additive_pq=pq_H,
        dense_H_eigenvalues=dense_eigenvalues,
        eta_nu_H_eigenvalues=structured_eigenvalues,
        additive_pq_H_eigenvalues=pq_eigenvalues,
    )

    summary = {
        "input": str(input_path.resolve()),
        "termination": metadata["termination"],
        "physical_covariance_dimension": int(covariance.shape[0]),
        "enlarged_frequency_dimension": int(frequency_factor.shape[0]),
        "eta_nu_real_latent_dimension": int(frequency_factor.shape[1]),
        "matsubara_takagi_tolerance": float(matsubara_tolerance),
        "physical_pseudocovariance_reconstruction_error": pseudo_error,
        "additive_pq_pseudocovariance_reconstruction_error": pq_pseudo_error,
        "eta_nu_zero_block_checks": {
            "M_nu_relative_maximum": float(max(x["M_nu_relative"] for x in block_data)),
            "nu_nu_relative_maximum": float(max(x["nu_nu_relative"] for x in block_data)),
        },
        "response_block_condition_number": {
            "minimum": float(min(x["condition"] for x in block_data)),
            "median": float(np.median([x["condition"] for x in block_data])),
            "maximum": float(max(x["condition"] for x in block_data)),
        },
        "Hermitian_covariance_comparison": {
            "trace_dense": dense_trace,
            "trace_eta_nu": structured_trace,
            "trace_ratio_eta_nu_over_dense": structured_trace / dense_trace,
            "largest_eigenvalue_dense": dense_max,
            "largest_eigenvalue_eta_nu": structured_max,
            "largest_eigenvalue_ratio_eta_nu_over_dense": structured_max / dense_max,
            "trace_additive_pq": pq_trace,
            "trace_ratio_additive_pq_over_dense": pq_trace / dense_trace,
            "largest_eigenvalue_additive_pq": pq_max,
            "largest_eigenvalue_ratio_additive_pq_over_dense": pq_max / dense_max,
            "relative_Frobenius_difference_additive_pq": float(pq_difference),
            "diagonal_ratio_additive_pq_minimum": float(np.min(pq_diagonal_ratio)),
            "diagonal_ratio_additive_pq_median": float(np.median(pq_diagonal_ratio)),
            "diagonal_ratio_additive_pq_maximum": float(np.max(pq_diagonal_ratio)),
            "relative_Frobenius_difference": float(difference),
            "diagonal_ratio_minimum": float(np.min(diagonal_ratio)),
            "diagonal_ratio_median": float(np.median(diagonal_ratio)),
            "diagonal_ratio_maximum": float(np.max(diagonal_ratio)),
            "effective_rank_dense": float(dense_trace**2 / np.sum(dense_eigenvalues**2)),
            "effective_rank_eta_nu": float(
                structured_trace**2 / np.sum(structured_eigenvalues**2)
            ),
            "modes_for_trace_fraction": {
                str(fraction): {
                    "dense": retained_rank(dense_eigenvalues, fraction),
                    "eta_nu": retained_rank(structured_eigenvalues, fraction),
                }
                for fraction in (0.9, 0.99, 0.999)
            },
        },
    }
    (output_dir / "keldysh_eta_nu_factor_comparison_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    analyze(args.input, args.output_dir)


if __name__ == "__main__":
    main()
