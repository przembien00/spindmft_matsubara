#!/usr/bin/env python3
"""Compare uncoupled spinDMFT (N=4, stype=A) vs Chebyshev (AFM rescale=0.5,
FM rescale=-0.5) for beta=0.5,1,1.5,2,2.5. One figure per correlation pair,
rows=beta, each panel shows all three datasets."""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import os

ROOT  = "/Users/przembien/Projects/spindmft_matsubara/CspinDMFT_finite_T"
DATA  = os.path.join(ROOT, "Data")
CHEB  = os.path.join(DATA, "Chebyshev")
PLOTS = os.path.join(ROOT, "Plots")

BETAS = [0.5, 1.0, 1.5, 2.0, 2.5]

# ── DMFT file map: beta → path (use first/positive run) ──────────────────────
def dmft_path(beta):
    b = str(beta).rstrip("0").rstrip(".")  if "." in str(beta) else str(beta)
    # beta as it appears in the filename (e.g. 0.5, 1, 1.5, 2, 2.5)
    bstr = f"{beta:g}"
    base = f"spinmodel=ISO__config=Square_2D_N=4_NN_J=0.0_uncoupled__uncoupled__beta={bstr}"
    # prefer the first (no suffix) file; for beta=2 that's the stype=B run –
    # use X (first stype=A run) instead
    for suffix in ["", "X", "XX"]:
        p = os.path.join(DATA, base + suffix + ".hdf5")
        if os.path.exists(p):
            # verify stype=A: shape[0]==1
            with h5py.File(p,"r") as f:
                if f["results/Re_correlation/0-0"].shape[0] == 1:
                    return p
    return None

def cheb_path(beta, rescale):
    bstr = f"{beta:g}"
    return os.path.join(CHEB, f"ISO__Square_NN_PBC_N=24__sites=0-1-7__beta={bstr}__rescale={rescale}.hdf5")

# ── loaders ───────────────────────────────────────────────────────────────────
def load_dmft(path, key):
    with h5py.File(path, "r") as f:
        ntp  = f["results/Re_correlation/0-0"].shape[1]
        tau  = np.linspace(0, 1, ntp)
        data = f["results/Re_correlation"][key][()][0]
        std  = f["runtimedata/Re_correlation_sample_stds"][key][()][0]
    return tau, data, std

def load_cheb(path, cheb_key):
    with h5py.File(path, "r") as f:
        ntp  = f["results/Re_correlation/0-0"].shape[1]
        data = f["results/Re_correlation"][cheb_key][()][0]
        std  = f["results/Re_stds"][cheb_key][()][0]
    # mirror around last datapoint to extend [0,beta/2] → [0,beta]
    data = np.concatenate([data, data[-2::-1]])
    std  = np.concatenate([std,  std[-2::-1]])
    tau  = np.linspace(0, 1, len(data))
    return tau, data, std

# ── pairs to plot ─────────────────────────────────────────────────────────────
cheb_key_map = {"0-0": "0-0", "0-1": "1-0", "0-3": "7-0"}
pairs = [
    ("0-0", r"$g^{xx}_{0-0}(\tau)$", "00"),
    ("0-1", r"$g^{xx}_{0-1}(\tau)$", "01"),
    ("0-3", r"$g^{xx}_{0-3}(\tau)$", "03"),
]

# ── style ─────────────────────────────────────────────────────────────────────
# Color by dataset type, marker by beta
dataset_colors = {
    "spinDMFT uncoupled": "#2ecc71",
    "Chebyshev AFM":      "#c0392b",
    "Chebyshev FM":       "#2980b9",
}
dataset_ls = {
    "spinDMFT uncoupled": "-",
    "Chebyshev AFM":      "--",
    "Chebyshev FM":       ":",
}
beta_markers = {0.5: "o", 1.0: "s", 1.5: "^", 2.0: "D", 2.5: "v"}

datasets = [
    ("spinDMFT uncoupled", "dmft", None),
    ("Chebyshev AFM",      "cheb", "0.5"),
    ("Chebyshev FM",       "cheb", "-0.5"),
]
MARKEVERY = 8

# ── one figure per correlation ────────────────────────────────────────────────
for dmft_key, ylabel, tag in pairs:
    cheb_key = cheb_key_map[dmft_key]

    fig, ax = plt.subplots(figsize=(7, 5))

    for beta in BETAS:
        mk = beta_markers[beta]
        for label, dtype, rescale in datasets:
            if dtype == "dmft":
                path = dmft_path(beta)
                if path is None:
                    continue
                tau, data, std = load_dmft(path, dmft_key)
            else:
                path = cheb_path(beta, rescale)
                if not os.path.exists(path):
                    continue
                tau, data, std = load_cheb(path, cheb_key)

            color = dataset_colors[label]
            ax.plot(tau, data, color=color, lw=1.3,
                    ls=dataset_ls[label], marker=mk, ms=4.5,
                    markevery=MARKEVERY, markerfacecolor="none")
            ax.fill_between(tau, data - std, data + std,
                            color=color, alpha=0.12)

    ax.axhline(0, color="gray", lw=0.6, ls=":")
    ax.set_xlim(0, 1)
    ax.set_xlabel(r"$\tau/\beta$", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)

    # legend: dataset type (colored lines) + beta (black markers only)
    import matplotlib.lines as mlines
    type_handles = [
        mlines.Line2D([], [], color=dataset_colors[lbl], ls=dataset_ls[lbl],
                      lw=1.5, label=lbl)
        for lbl in dataset_colors
    ]
    beta_handles = [
        mlines.Line2D([], [], color="k", ls="none",
                      marker=beta_markers[b], ms=5, markerfacecolor="none",
                      label=fr"$\beta={b:g}$")
        for b in BETAS
    ]
    ax.legend(handles=type_handles + beta_handles,
              fontsize=9, framealpha=0.9, loc="best", ncol=1)

    plt.tight_layout()
    out_base = os.path.join(PLOTS, f"uncoupled_vs_chebyshev_N4_{tag}")
    for ext in ("pdf", "png"):
        plt.savefig(f"{out_base}.{ext}", dpi=150 if ext == "png" else None,
                    bbox_inches="tight")
        print(f"Saved: {out_base}.{ext}")
    plt.close()

print("Done.")
