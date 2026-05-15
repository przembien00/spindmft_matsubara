#!/usr/bin/env python3
"""Compare N=4 beta=1 correlations: spinDMFT (coupled J=±0.5 & uncoupled)
vs Chebyshev benchmark. One file per correlation: 0-0, 0-1, 0-3."""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import os

# ── paths ────────────────────────────────────────────────────────────────────
ROOT  = "/Users/przembien/Projects/spindmft_matsubara/CspinDMFT_finite_T"
DATA  = os.path.join(ROOT, "Data")
CHEB  = os.path.join(DATA, "Chebyshev")
BETA  = 1

# ── DMFT files ───────────────────────────────────────────────────────────────
dmft_files = {
    "J=0.5"           : os.path.join(DATA, "spinmodel=ISO__config=Square_2D_N=4_NN_J=0.5__beta=1.hdf5"),
    "J=-0.5"          : os.path.join(DATA, "spinmodel=ISO__config=Square_2D_N=4_NN_J=-0.5__beta=1.hdf5"),
    "J=0 (uncoupled)" : os.path.join(DATA, "spinmodel=ISO__config=Square_2D_N=4_NN_J=0.0_uncoupled__uncoupled__beta=1.hdf5"),
}

# ── Chebyshev files ───────────────────────────────────────────────────────────
cheb_files = {
    "J=0.5"  : os.path.join(CHEB, f"ISO__Square_NN_PBC_N=24__sites=0-1-7__beta={BETA}__rescale=0.5.hdf5"),
    "J=-0.5" : os.path.join(CHEB, f"ISO__Square_NN_PBC_N=24__sites=0-1-7__beta={BETA}__rescale=-0.5.hdf5"),
}

# Chebyshev dataset key → DMFT pair label
cheb_key_map = {"0-0": "0-0", "1-0": "0-1", "7-0": "0-3"}
dmft_to_cheb = {v: k for k, v in cheb_key_map.items()}

# ── style ─────────────────────────────────────────────────────────────────────
# Colors: one per J value (shared between DMFT and Chebyshev for easy comparison)
colors = {
    "J=0.5"           : "#c0392b",   # red  – AFM
    "J=-0.5"          : "#2980b9",   # blue – FM
    "J=0 (uncoupled)" : "#27ae60",   # green
}
# Linestyle and marker distinguish method
dmft_style  = dict(ls="-",  marker="o", ms=4,  lw=1.5, markevery=5,  markerfacecolor="none")
cheb_style  = dict(ls="--", marker="s", ms=4,  lw=1.5, markevery=10, markerfacecolor="none")

pairs = [
    ("0-0", r"on-site $g^{ab}(\tau)$",  "00"),
    ("0-1", r"NN $g^{ab}(\tau)$",        "01"),
    ("0-3", r"NNN $g^{ab}(\tau)$",       "03"),
]

# ── load helpers ──────────────────────────────────────────────────────────────
def load_dmft(path, key):
    with h5py.File(path, "r") as f:
        ntp  = f["results/Re_correlation/0-0"].shape[1]
        tau  = np.linspace(0, 1, ntp)
        data = f["results/Re_correlation"][key][()][0]
        std  = f["runtimedata/Re_correlation_sample_stds"][key][()][0]
    return tau, data, std

def load_cheb(path, key):
    """Load Chebyshev data (defined on [0, beta/2]) and mirror around the
    last datapoint to extend to [0, beta]."""
    with h5py.File(path, "r") as f:
        ntp  = f["results/Re_correlation/0-0"].shape[1]   # points on [0, beta/2]
        data = f["results/Re_correlation"][key][()][0]
        std  = f["results/Re_stds"][key][()][0]
    # mirror: append reversed data (excluding the last point) after the last point
    data_full = np.concatenate([data, data[-2::-1]])
    std_full  = np.concatenate([std,  std[-2::-1]])
    tau_full  = np.linspace(0, 1, len(data_full))
    return tau_full, data_full, std_full

# ── one figure per correlation ────────────────────────────────────────────────
for pair, ylabel, tag in pairs:
    fig, ax = plt.subplots(figsize=(6, 4.5))

    # DMFT
    for label, path in dmft_files.items():
        if not os.path.exists(path):
            continue
        tau, data, std = load_dmft(path, pair)
        ax.plot(tau, data, color=colors[label],
                label=f"spinDMFT {label}", **dmft_style)
        ax.fill_between(tau, data - std, data + std,
                        color=colors[label], alpha=0.15)

    # Chebyshev
    cheb_key = dmft_to_cheb[pair]
    for label, path in cheb_files.items():
        if not os.path.exists(path):
            continue
        tau, data, std = load_cheb(path, cheb_key)
        ax.plot(tau, data, color=colors[label],
                label=f"Chebyshev {label}", **cheb_style)
        ax.fill_between(tau, data - std, data + std,
                        color=colors[label], alpha=0.15)

    ax.axhline(0, color="gray", lw=0.7, ls=":")
    ax.set_xlim(0, 1)
    if pair == "0-0":
        ax.set_ylim(0.2, 0.25)
    ax.set_xlabel(r"$\tau/\beta$", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.legend(fontsize=9, ncol=1, framealpha=0.9)

    plt.tight_layout()
    out_base = os.path.join(ROOT, f"plot_comparison_N4_beta{BETA}_{tag}")
    for ext in ("pdf", "png"):
        plt.savefig(f"{out_base}.{ext}", dpi=150 if ext == "png" else None,
                    bbox_inches="tight")
        print(f"Saved: {out_base}.{ext}")
    plt.close()
