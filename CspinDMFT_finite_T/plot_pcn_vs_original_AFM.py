#!/usr/bin/env python3
"""Compare the pCN sampler (Algorithm_MH) against the original brute-force
algorithm for the N=2 AFM cluster (ISO, J=0.5) at beta = 1,2,3,4,5.

All curves are drawn on a single axis: colour encodes beta, line style encodes
the method (solid = original importance sampling, dashed = pCN). The on-site
imaginary-time correlation g^{aa}_{00}(tau) is shown, with tau normalised by beta.
"""

import os
import glob
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

ROOT = "/Users/przembien/Projects/spindmft_matsubara/CspinDMFT_finite_T"
matplotlib.rc_file(os.path.join(ROOT, "matplotlibrc"))

DATA  = os.path.join(ROOT, "Data")
PLOTS = os.path.join(ROOT, "Plots")
os.makedirs(PLOTS, exist_ok=True)

BETAS = [1, 2, 3, 4, 5]
PAIR  = "0-0"   # on-site correlation

def original_path(beta):
    return os.path.join(DATA, "AFM",
        f"spinmodel=ISO__config=Square_2D_N=2_NN_J=0.5__beta={beta}.hdf5")

def pcn_path(beta):
    hits = glob.glob(os.path.join(DATA, "**",
        f"spinmodel=ISO__config=Square_2D_N=2_NN_J=0.5__beta={beta}_AFM_pcn.hdf5"),
        recursive=True)
    return hits[0] if hits else None

def load(path, pair):
    with h5py.File(path, "r") as f:
        data = f["results/Re_correlation"][pair][()][0]   # stype A: single direction
        tau  = np.linspace(0.0, 1.0, data.shape[0])       # tau / beta in [0, 1]
    return tau, data

# one colour per beta
cmap   = plt.cm.viridis
colors = {b: cmap(i / (len(BETAS) - 1)) for i, b in enumerate(BETAS)}

fig, ax = plt.subplots(figsize=(6.0, 4.2))

for b in BETAS:
    tau_o, g_o = load(original_path(b), PAIR)
    ax.plot(tau_o, g_o, ls="-", color=colors[b], lw=1.5,
            label=rf"$\beta={b}$")

    p = pcn_path(b)
    if p is not None:
        tau_p, g_p = load(p, PAIR)
        ax.plot(tau_p, g_p, ls="--", color=colors[b], lw=1.5,
                marker="o", ms=4, markevery=8, markerfacecolor="none")

# method legend (line-style), kept separate from the beta (colour) legend
orig_handle = plt.Line2D([], [], color="black", ls="-",  lw=1.5, label="original")
pcn_handle  = plt.Line2D([], [], color="black", ls="--", lw=1.5,
                         marker="o", ms=4, markerfacecolor="none", label="pCN")

leg_beta = ax.legend(loc="upper center", ncol=5, columnspacing=1.0,
                     handlelength=1.4, title=None)
ax.add_artist(leg_beta)
ax.legend(handles=[orig_handle, pcn_handle], loc="lower center")

ax.set_xlim(0.0, 1.0)
ax.set_xlabel(r"$\tau / \beta$")
ax.set_ylabel(r"$g^{aa}_{00}(\tau)$")

out = os.path.join(PLOTS, "pcn_vs_original_AFM_N2_g00.pdf")
fig.savefig(out)
fig.savefig(out.replace(".pdf", ".png"), dpi=200)
print("saved:", out)
