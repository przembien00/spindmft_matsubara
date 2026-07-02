#!/usr/bin/env python3
"""Before/after validation of the pCN observable-contribution cache (rejected
steps re-accumulate cached values instead of recomputing the propagator
products). Both builds ran the N=2 cluster (ISO, J=0.5, stype=A) at
beta = 1,2,3 with the same seed, so the curves must coincide.

Colour encodes beta, line style encodes the build (solid = baseline,
dashed + open markers = cached). Error bars are the blocking standard errors.
Top panel: on-site g^{aa}_{00}(tau); bottom panel: inter-site g^{aa}_{01}(tau);
tau is normalised by beta.
"""

import os
import glob
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.abspath(__file__))
matplotlib.rc_file(os.path.join(ROOT, "matplotlibrc"))

DATA  = os.path.join(ROOT, "Data", "SpeedupBench")
PLOTS = os.path.join(ROOT, "Plots")
os.makedirs(PLOTS, exist_ok=True)

BETAS = [1, 2, 3]
PAIRS = ["0-0", "0-1"]


def find(beta, build):
    pattern = os.path.join(
        DATA, "**",
        f"spinmodel=ISO__config=Square_2D_N=2_NN_J=0.5__beta={beta}_{build}_b{beta}.hdf5")
    hits = glob.glob(pattern, recursive=True)
    if not hits:
        raise FileNotFoundError(pattern)
    return hits[0]


def load(path, pair):
    with h5py.File(path, "r") as f:
        g   = f["results/Re_correlation"][pair][()][0]                  # stype A: single direction
        std = f["runtimedata/Re_correlation_sample_stds"][pair][()][0]
        tau = np.linspace(0.0, 1.0, g.shape[0])                          # tau / beta in [0, 1]
    return tau, g, std


cmap   = plt.cm.viridis
colors = {b: cmap(i / (len(BETAS) - 1)) for i, b in enumerate(BETAS)}

fig, axes = plt.subplots(2, 1, figsize=(6.0, 6.6), sharex=True)

for ax, pair in zip(axes, PAIRS):
    for b in BETAS:
        tau, g, std = load(find(b, "baseline"), pair)
        ax.errorbar(tau, g, yerr=std, ls="-", color=colors[b], lw=1.5,
                    errorevery=5, capsize=2, label=rf"$\beta={b}$")

        tau, g, std = load(find(b, "cached"), pair)
        ax.errorbar(tau, g, yerr=std, ls="--", color=colors[b], lw=1.5,
                    marker="o", ms=4, markevery=7, markerfacecolor="none",
                    errorevery=(3, 5), capsize=2)
    ax.set_xlim(0.0, 1.0)

base_handle   = plt.Line2D([], [], color="black", ls="-",  lw=1.5, label="baseline")
cached_handle = plt.Line2D([], [], color="black", ls="--", lw=1.5,
                           marker="o", ms=4, markerfacecolor="none", label="cached")

leg_beta = axes[0].legend(loc="upper center", ncol=3, columnspacing=1.0,
                          handlelength=1.4)
axes[0].add_artist(leg_beta)
axes[0].legend(handles=[base_handle, cached_handle], loc="lower center")

axes[0].set_ylabel(r"$g^{aa}_{00}(\tau)$")
axes[1].set_ylabel(r"$g^{aa}_{01}(\tau)$")
axes[1].set_xlabel(r"$\tau / \beta$")

out = os.path.join(PLOTS, "cache_speedup_comparison_N2.pdf")
fig.savefig(out)
fig.savefig(out.replace(".pdf", ".png"), dpi=200)
print("saved:", out)

# numeric confirmation of the max deviation between the two builds
for pair in PAIRS:
    for b in BETAS:
        _, g0, s0 = load(find(b, "baseline"), pair)
        _, g1, s1 = load(find(b, "cached"),   pair)
        print(f"pair {pair} beta={b}: max|dg|={np.max(np.abs(g1-g0)):.3e}  "
              f"max|dstd|={np.max(np.abs(s1-s0)):.3e}")
