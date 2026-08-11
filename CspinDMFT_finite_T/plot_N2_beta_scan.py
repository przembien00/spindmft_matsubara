#!/usr/bin/env python3
"""
Plot imaginary-time spin correlations from N=2 uncoupled spins beta scan.
Creates separate figures for site pairs (0,0), (0,1),
each overlaying all 7 beta values (1, 2, 3, 5, 8, 12, 18).

Correlation data: Re_correlation (real part of g^xx_{ij}(tau))
Color: colormap by beta
"""

import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Load matplotlibrc from project directory
plt.style.use(os.path.join(os.path.dirname(__file__), "matplotlibrc"))

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(ROOT, "Data", "Uncoupled_N2_scan")
PLOTS = os.path.join(ROOT, "Plots", "Uncoupled_N2_beta_scan")
os.makedirs(PLOTS, exist_ok=True)

# Betas to plot
BETAS = [1, 2, 3, 5, 8, 12, 18]

# Site pairs to create separate figures for
PAIRS = [
    ("0-0", r"$g^{xx}_{0,0}(\tau)$"),
    ("0-1", r"$g^{xx}_{0,1}(\tau)$"),
]

# Color mapping by beta
CMAP = plt.get_cmap("plasma_r")
BETA_NORM = mcolors.Normalize(vmin=min(BETAS), vmax=max(BETAS))

def beta_color(beta):
    """Return color for a given beta value."""
    return CMAP(BETA_NORM(beta))

def data_path(beta):
    """Construct path to data file for a given beta."""
    bstr = f"{beta:g}"
    # Check for both base filename and _ETV variant
    base = f"spinmodel=ISO__config=Square_2D_N=2_NN_J=0.0_uncoupled__beta={bstr}"
    for suffix in ["", "_ETV"]:
        path = os.path.join(DATA, base + suffix + ".hdf5")
        if os.path.exists(path):
            return path
    return None

def load_correlation(path, pair_key):
    """
    Load correlation data from HDF5 file.

    Returns:
        tau: imaginary-time grid normalized to [0, 1]
        data: correlation values
        stds: statistical error bars
    """
    with h5py.File(path, "r") as f:
        ntp = f["results/Re_correlation"][pair_key].shape[1]

        # Load real part (physical correlations)
        data = f["results/Re_correlation"][pair_key][()][0].astype(float)
        stds = f["runtimedata/Re_correlation_sample_stds"][pair_key][()][0].astype(float)

    # Normalize tau to [0, 1]
    tau = np.linspace(0, 1, ntp)

    return tau, data, stds

# Create one figure per site pair
for pair_key, ylabel in PAIRS:
    fig, ax = plt.subplots(figsize=(7, 5))

    ymin, ymax = np.inf, -np.inf

    # Plot each beta
    for beta in BETAS:
        path = data_path(beta)
        if path is None:
            print(f"Warning: data file not found for beta={beta}")
            continue

        try:
            tau, data, stds = load_correlation(path, pair_key)
        except (KeyError, OSError) as e:
            print(f"Warning: could not load {pair_key} from beta={beta}: {e}")
            continue

        color = beta_color(beta)

        # Plot correlation with error band
        ax.plot(tau, data, color=color, lw=1.3, label=fr"$\beta={beta}$")
        ax.fill_between(tau, data - stds, data + stds, color=color, alpha=0.15)

        ymin = min(ymin, np.min(data - stds))
        ymax = max(ymax, np.max(data + stds))

    # Style axes
    ax.axhline(0, color="gray", lw=0.6, ls=":")
    ax.set_xlim(0, 1)
    ax.margins(x=0)  # Remove x-axis padding

    if np.isfinite(ymin) and np.isfinite(ymax):
        margin = 0.05 * (ymax - ymin)
        ax.set_ylim(ymin - margin, ymax + margin)

    ax.set_xlabel(r"$\tau/\beta$", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)

    # Legend: beta values with colors
    ax.legend(fontsize=10, loc="best", framealpha=0.9)

    plt.tight_layout()

    # Save with descriptive name based on site pair
    pair_tag = pair_key.replace("-", "")
    out_file = os.path.join(PLOTS, f"g_xx_{pair_tag}.pdf")
    plt.savefig(out_file, bbox_inches="tight")
    print(f"Saved: {out_file}")

    plt.close()

print("Done.")
