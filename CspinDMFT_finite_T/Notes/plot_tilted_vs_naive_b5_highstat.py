#!/usr/bin/env python3
"""High-statistics check at beta=5: naive vs Laplace-tilted IS for N=2, J=0.5.
4 cores x 100000 samples = 400000 samples per iteration."""

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

PAIRS = [("0-0", r"$g^{xx}_{ii}(\tau)$"),
         ("0-1", r"$g^{xx}_{ij}(\tau)$")]

def load(tag):
    fn = os.path.join(DATA,
        f"spinmodel=ISO__config=Square_2D_N=2_NN_J=0.5__beta=5_{tag}.hdf5")
    with h5py.File(fn, "r") as f:
        out = {}
        for k, _ in PAIRS:
            corr = f["results/Re_correlation"][k][()][0]
            std  = f["runtimedata/Re_correlation_sample_stds"][k][()][0]
            out[k] = (corr, std)
        ess = np.asarray(f["runtimedata"].attrs.get("ess_relative_list", np.array([])))
    return out, ess

naive,  ess_n = load("ESS_b5_M400k_naive")
tilted, ess_t = load("ESS_b5_M400k_tilted")

fig, axes = plt.subplots(1, len(PAIRS), figsize=(9.5, 4.0))

colors = {"naive": "#c0392b", "tilted": "#2980b9"}
labels = {"naive":  "naive $q=p$",
          "tilted": r"Laplace-tilted $q\!\propto\!p\cdot Z_{\rm Laplace}$"}

tau = np.linspace(0, 1, naive["0-0"][0].size)
step = 4

for col, (key, ylabel) in enumerate(PAIRS):
    ax = axes[col]
    for offset, (tag, dd) in enumerate((("naive", naive), ("tilted", tilted))):
        data, std = dd[key]
        c = colors[tag]
        idx = np.arange(offset, tau.size, step)
        ax.plot(tau, data, color=c, lw=1.2, label=labels[tag])
        ax.errorbar(tau[idx], data[idx], yerr=std[idx],
                    fmt="none", ecolor=c, elinewidth=0.9,
                    capsize=2.2, capthick=0.9, alpha=0.9)
    ax.set_xlim(0, 1)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(r"$\tau/\beta$")

ess_n_plat = float(np.mean(ess_n[-3:])) if ess_n.size else float("nan")
ess_t_plat = float(np.mean(ess_t[-3:])) if ess_t.size else float("nan")
axes[0].text(0.03, 0.05,
    rf"$\beta=5$, $M=400{{,}}000$"
    "\n"
    rf"ESS/M: naive {ess_n_plat:.2f},  tilted {ess_t_plat:.2f}",
    transform=axes[0].transAxes, va="bottom", ha="left",
    fontsize=9,
    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=2))
axes[0].legend(loc="upper right", fontsize=9, framealpha=0.9)

plt.tight_layout()
for ext in ("pdf", "png"):
    p = os.path.join(OUT, f"tilted_vs_naive_b5_M400k.{ext}")
    plt.savefig(p, dpi=160 if ext == "png" else None, bbox_inches="tight")
    print(f"Saved: {p}")
plt.close()

# Also print a numerical diff diagnostic
print("\nMean-curve discrepancy (tilted - naive) in units of combined std:")
for key, _ in PAIRS:
    d_n, s_n = naive[key]; d_t, s_t = tilted[key]
    diff = d_t - d_n
    sigma = np.sqrt(s_n**2 + s_t**2)
    z = diff / sigma
    print(f"  {key}: |z|_max = {np.max(np.abs(z)):.2f},  mean |z| = {np.mean(np.abs(z)):.2f}")
