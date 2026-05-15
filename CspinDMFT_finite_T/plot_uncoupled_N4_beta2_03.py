#!/usr/bin/env python3
"""Compare 0-3 (NNN) correlation for N=4 beta=2 uncoupled with
initpaircorr=positive vs negative (stype=A)."""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import os

DATA = "/Users/przembien/Projects/spindmft_matsubara/CspinDMFT_finite_T/Data"
BETA = 2.0

files = {
    "positive": os.path.join(DATA, "spinmodel=ISO__config=Square_2D_N=4_NN_J=0.0_uncoupled__uncoupled__beta=2X.hdf5"),
    "negative": os.path.join(DATA, "spinmodel=ISO__config=Square_2D_N=4_NN_J=0.0_uncoupled__uncoupled__beta=2XX.hdf5"),
}
colors = {"positive": "steelblue", "negative": "firebrick"}

fig, ax = plt.subplots(figsize=(6, 4.5))

for label, path in files.items():
    with h5py.File(path, "r") as f:
        ntp  = f["results/Re_correlation/0-0"].shape[1]
        tau  = np.linspace(0, 1, ntp)
        data = f["results/Re_correlation"]["0-3"][()][0]
        std  = f["runtimedata/Re_correlation_sample_stds"]["0-3"][()][0]
    ax.plot(tau, data, ls="-", marker="o", ms=3, color=colors[label],
            label=f"init={label}")
    ax.fill_between(tau, data - std, data + std, color=colors[label], alpha=0.2)

ax.axhline(0, color="gray", lw=0.7, ls=":")
ax.set_xlim(0, 1)
ax.set_xlabel(r"$\tau/\beta$", fontsize=12)
ax.set_ylabel(r"$g^{ab}_{03}(\tau)$", fontsize=12)
ax.legend(fontsize=10)

plt.tight_layout()
ROOT = "/Users/przembien/Projects/spindmft_matsubara/CspinDMFT_finite_T"
for ext in ("pdf", "png"):
    out = os.path.join(ROOT, f"plot_uncoupled_N4_beta2_03.{ext}")
    plt.savefig(out, dpi=150 if ext == "png" else None, bbox_inches="tight")
    print(f"Saved: {out}")
plt.show()
