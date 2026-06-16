#!/usr/bin/env python3
"""Compare naive vs Laplace-tilted importance sampling for N=2, J=0.5
across beta = 1, 3, 5.  Plots Re g^{ab}_{ij}(tau) with sample std bands."""

import os
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = "/Users/przembien/Projects/spindmft_matsubara/CspinDMFT_finite_T"
mpl.rc_file(os.path.join(ROOT, "matplotlibrc"))

DATA = os.path.join(ROOT, "Data", "ESS_diagnostic")
OUT  = os.path.join(ROOT, "Plots")
os.makedirs(OUT, exist_ok=True)

BETAS = [1, 3, 5]
PAIRS = [("0-0", r"$g^{xx}_{ii}(\tau)$"),
         ("0-1", r"$g^{xx}_{ij}(\tau)$")]

def load(beta, tag):
    fn = os.path.join(
        DATA,
        f"spinmodel=ISO__config=Square_2D_N=2_NN_J=0.5__beta={beta}_{tag}.hdf5",
    )
    with h5py.File(fn, "r") as f:
        out = {}
        for k, _ in PAIRS:
            corr = f["results/Re_correlation"][k][()][0]  # xx channel
            std  = f["runtimedata/Re_correlation_sample_stds"][k][()][0]
            out[k] = (corr, std)
        ess = f["runtimedata"].attrs.get("ess_relative_list", np.array([]))
    return out, np.asarray(ess)

fig, axes = plt.subplots(len(BETAS), len(PAIRS), figsize=(8.5, 8.0),
                        sharex=True)

colors = {"naive": "#c0392b", "tilted": "#2980b9"}
labels = {"naive":  "naive $q=p$",
          "tilted": r"Laplace-tilted $q\!\propto\!p\cdot Z_{\rm Laplace}$"}

for row, beta in enumerate(BETAS):
    naive,  ess_n = load(beta, "ESS_full")
    tilted, ess_t = load(beta, "ESS_tilted_full")
    tau = np.linspace(0, 1, naive["0-0"][0].size)

    for col, (key, ylabel) in enumerate(PAIRS):
        ax = axes[row, col]
        for offset, (tag, datadict) in enumerate((("naive", naive), ("tilted", tilted))):
            data, std = datadict[key]
            c = colors[tag]
            # downsample errorbars so they remain legible
            step = 6
            idx = np.arange(offset, tau.size, step)
            ax.plot(tau, data, color=c, lw=1.2, label=labels[tag])
            ax.errorbar(tau[idx], data[idx], yerr=std[idx],
                        fmt="none", ecolor=c, elinewidth=0.9,
                        capsize=2.0, capthick=0.9, alpha=0.9)
        ax.set_xlim(0, 1)
        if col == 0:
            ax.set_ylabel(ylabel)
        if row == len(BETAS) - 1:
            ax.set_xlabel(r"$\tau/\beta$")
        # Annotate beta and converged ESS/M
        ess_n_plat = float(np.mean(ess_n[-3:])) if ess_n.size else float("nan")
        ess_t_plat = float(np.mean(ess_t[-3:])) if ess_t.size else float("nan")
        ax.text(0.03, 0.95,
                rf"$\beta={beta}$" "\n"
                rf"ESS/M: naive {ess_n_plat:.2f},  tilted {ess_t_plat:.2f}",
                transform=ax.transAxes, va="top", ha="left",
                fontsize=9,
                bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=2))

# one legend at top
axes[0, 0].legend(loc="lower right", fontsize=9, framealpha=0.9)

plt.tight_layout()
for ext in ("pdf", "png"):
    out_path = os.path.join(OUT, f"tilted_vs_naive_N2_J0.5.{ext}")
    plt.savefig(out_path, dpi=160 if ext == "png" else None,
                bbox_inches="tight")
    print(f"Saved: {out_path}")
plt.close()
