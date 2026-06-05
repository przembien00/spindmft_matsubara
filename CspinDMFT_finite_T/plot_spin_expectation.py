#!/usr/bin/env python3
"""Quiver plots of thermal spin expectation values on the 3x3 lattice.

One figure per magnetic type (AFM / FM), with one panel per beta.
Arrows show (Sx, Sy); color encodes Sz. Shared colorbar per figure.
"""

import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection

ROOT  = os.path.dirname(os.path.abspath(__file__))
PLOTS = os.path.join(ROOT, "Plots", "spin_expectation")
os.makedirs(PLOTS, exist_ok=True)

BETAS = [1, 2, 3, 4, 5]
CMAP  = plt.get_cmap("RdBu_r")

# grid dimensions per N
GRID = {2: (2, 1), 4: (2, 2), 9: (3, 3)}

def make_lattice(Lx, Ly):
    N   = Lx * Ly
    pos = np.array([(i % Lx, i // Lx) for i in range(N)], dtype=float)
    bonds = [
        (pos[i], pos[j])
        for i in range(N) for j in range(i + 1, N)
        if abs(pos[i, 0] - pos[j, 0]) + abs(pos[i, 1] - pos[j, 1]) == 1
    ]
    return pos, bonds

def data_path(mag_type, N, beta):
    J = 0.5 if mag_type == "AFM" else -0.5
    return os.path.join(ROOT, "Data", mag_type,
        f"spinmodel=ISO__config=Square_2D_N={N}_NN_J={J}__beta={beta:g}.hdf5")

for N, (Lx, Ly) in GRID.items():
    pos, bonds = make_lattice(Lx, Ly)
    # panel aspect ratio: slightly taller than the grid to accommodate title
    panel_h = max(Ly / Lx * 3.0, 1.8)

    for mag_type in ("AFM", "FM"):
        # global Sz range for colorbar reference
        all_sz = []
        for beta in BETAS:
            p = data_path(mag_type, N, beta)
            if os.path.exists(p):
                with h5py.File(p, "r") as f:
                    all_sz.append(f["results/spin_expectation_values"][()][:, 2])
        if not all_sz:
            continue
        vlim = np.max(np.abs(np.concatenate(all_sz))) * 1.05
        norm = mcolors.Normalize(vmin=-vlim, vmax=vlim)

        fig = plt.figure(figsize=(3 * len(BETAS) + 0.8, panel_h + 0.6))
        gs  = fig.add_gridspec(1, len(BETAS) + 1,
                               width_ratios=[1] * len(BETAS) + [0.08],
                               wspace=0.05)
        axes    = [fig.add_subplot(gs[0, k]) for k in range(len(BETAS))]
        cbar_ax = fig.add_subplot(gs[0, len(BETAS)])

        for ax, beta in zip(axes, BETAS):
            p = data_path(mag_type, N, beta)
            ax.set_title(rf"$\beta J_Q={beta:g}$", fontsize=10)

            if not os.path.exists(p):
                ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                        transform=ax.transAxes, fontsize=12, color="gray")
                ax.axis("off")
                continue

            with h5py.File(p, "r") as f:
                sev = f["results/spin_expectation_values"][()]
            Sx, Sy, Sz = sev[:, 0], sev[:, 1], sev[:, 2]

            xy_max     = np.hypot(Sx, Sy).max()
            scale      = 0.4 / xy_max if xy_max > 0 else 1.0
            panel_vlim = np.max(np.abs(Sz)) * 1.05
            panel_norm = mcolors.Normalize(vmin=-panel_vlim, vmax=panel_vlim)

            lc = LineCollection([np.array([p0, p1]) for p0, p1 in bonds],
                                colors="0.80", linewidths=1.2, zorder=1)
            ax.add_collection(lc)
            ax.scatter(pos[:, 0], pos[:, 1], c=Sz, cmap=CMAP, norm=panel_norm,
                       s=60, zorder=2, edgecolors="0.3", linewidths=0.5)
            ax.quiver(pos[:, 0], pos[:, 1], Sx * scale, Sy * scale,
                      color="0.15", angles="xy", scale_units="xy", scale=1,
                      width=0.008, headwidth=5, headlength=5, headaxislength=4,
                      zorder=3)

            ax.set_xlim(-0.6, Lx - 0.4)
            ax.set_ylim(-0.6, Ly - 0.4)
            ax.set_aspect("equal")
            ax.axis("off")

        sm = plt.cm.ScalarMappable(cmap=CMAP, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, cax=cbar_ax).set_label(
            r"$\langle S^z_i \rangle$", fontsize=10)

        fname = os.path.join(PLOTS, f"spin_expectation_quiver_{mag_type}_N{N}_betas.pdf")
        plt.savefig(fname, bbox_inches="tight")
        print(f"Saved {fname}")
        plt.close()
