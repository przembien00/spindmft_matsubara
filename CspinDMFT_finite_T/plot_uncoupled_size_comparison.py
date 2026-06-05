#!/usr/bin/env python3
"""Compare spinDMFT (N=2,4,9) and Chebyshev (N=24, rescale=±0.5) for the uncoupled (J=0) case.

Single figure with three vertically stacked axes (local, NN, NNN) sharing the same x axis.
Color = beta (plasma_r: yellow=small, purple=large).
Line style = N for DMFT; solid + triangle markers for Chebyshev.
Marker: ^ for rescale=+0.5, v for rescale=-0.5.
"""

import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
import matplotlib.colors as mcolors

plt.rcParams["figure.autolayout"] = False

ROOT      = os.path.dirname(os.path.abspath(__file__))
DATA_UNC  = os.path.join(ROOT, "Data", "Uncoupled")
DATA_CHEB = os.path.join(ROOT, "Data", "Chebyshev")
PLOTS     = os.path.join(ROOT, "Plots", "uncoupled_size_comparison")
os.makedirs(PLOTS, exist_ok=True)

BETAS = [0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0][::2]

N_LSTYLE       = {2: "--", 4: "-.", 9: ":"}
RESCALE_MARKER = {0.5: "^", -0.5: "v"}

CMAP      = plt.get_cmap("plasma_r")
BETA_NORM = mcolors.Normalize(vmin=min(BETAS), vmax=max(BETAS))

def beta_color(beta: float):
    return CMAP(BETA_NORM(beta))

def uncoupled_path(N: int, beta: float) -> str:
    bstr = f"{beta:g}"
    if N == 2:
        fname = f"spinmodel=ISO__config=Square_2D_N={N}_NN_J=0.0_uncoupled__beta={bstr}.hdf5"
    else:
        fname = f"spinmodel=ISO__config=Square_2D_N={N}_NN_J=0.0_uncoupled__uncoupled__beta={bstr}.hdf5"
    return os.path.join(DATA_UNC, fname)

def cheb_path(rescale: float, beta: float) -> str:
    bstr = f"{beta:g}"
    return os.path.join(DATA_CHEB,
        f"ISO__Square_NN_PBC_N=24__sites=0-1-7__beta={bstr}__rescale={rescale}.hdf5")

def load_dmft(path: str, key: str):
    with h5py.File(path, "r") as f:
        delta_t = float(f["parameters"].attrs["delta_t"])
        ntp     = f["results/Re_correlation/0-0"].shape[1]
        data    = f["results/Re_correlation"][key][()][0].astype(float)
        stds    = f["runtimedata/Re_correlation_sample_stds"][key][()][0].astype(float)
        n_samp  = float(f["parameters"].attrs["num_Samples"])
    tau = np.arange(ntp) * delta_t / (delta_t * (ntp - 1))
    return tau, data, stds / np.sqrt(n_samp)

def load_cheb(path: str, key: str):
    with h5py.File(path, "r") as f:
        data = f["results/Re_correlation"][key][()][0].astype(float)
        std  = f["results/Re_stds"][key][()][0].astype(float)
    data = np.concatenate([data, data[-2::-1]])
    std  = np.concatenate([std,  std[-2::-1]])
    tau  = np.linspace(0, 1, len(data))
    return tau, data, std

CORR_SPECS = [
    dict(
        name="local",
        dmft_keys={2: "0-0", 4: "0-0", 9: "0-0"},
        cheb_key="0-0",
        ylabel=r"$g^{xx}_\mathrm{local}(\tau)$",
    ),
    dict(
        name="NN",
        dmft_keys={2: "0-1", 4: "0-1", 9: "0-1"},
        cheb_key="1-0",
        ylabel=r"$g^{xx}_\mathrm{NN}(\tau)$",
    ),
    dict(
        name="NNN",
        dmft_keys={2: None, 4: "0-3", 9: "0-4"},
        cheb_key="7-0",
        ylabel=r"$g^{xx}_\mathrm{NNN}(\tau)$",
    ),
]

MARKER_EVERY = 10

fig, axes = plt.subplots(3, 1, figsize=(7, 7.5), sharex=True,
                         gridspec_kw={"hspace": 0})
fig.subplots_adjust(left=0.11, right=0.85, top=0.95, bottom=0.09)

for row, (ax, spec) in enumerate(zip(axes, CORR_SPECS)):
    ymin, ymax = np.inf, -np.inf

    for beta in BETAS:
        color = beta_color(beta)

        # ── DMFT N=2,4,9 ──
        for N in (2, 4, 9):
            dmft_key = spec["dmft_keys"][N]
            if dmft_key is None:
                continue
            path = uncoupled_path(N, beta)
            if not os.path.exists(path):
                continue
            try:
                tau, data, err = load_dmft(path, dmft_key)
            except (KeyError, OSError):
                continue
            ax.plot(tau, data, color=color, lw=1.2, ls=N_LSTYLE[N])
            ax.fill_between(tau, data - err, data + err,
                            color=color, alpha=0.12)
            ymin = min(ymin, np.min(data))
            ymax = max(ymax, np.max(data))


    ax.axhline(0, color="gray", lw=0.5, ls=":")
    ax.set_xlim(0, 1)
    if np.isfinite(ymin) and np.isfinite(ymax):
        margin = 0.05 * (ymax - ymin)
        ax.set_ylim(ymin - margin, ymax + margin)
    ax.set_ylabel(spec["ylabel"])

    if row < 2:
        ax.tick_params(labelbottom=False)

axes[-1].set_xlabel(r"$\tau/\beta$")

# Legend on top axis
n_handles = []
for N, ls in N_LSTYLE.items():
    n_handles.append(mlines.Line2D([], [], color="k", ls=ls,
                                   lw=1.5, label=f"$N={N}$ (spinDMFT)"))
axes[0].legend(handles=n_handles, fontsize=8, loc="lower left", framealpha=0.85)

# Shared colorbar placed to the right of all axes
cbar_ax = fig.add_axes([0.88, 0.09, 0.025, 0.95 - 0.09])
sm = cm.ScalarMappable(cmap=CMAP, norm=BETA_NORM)
sm.set_array([])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label(r"$\beta J_Q$", fontsize=13)
cbar.set_ticks(BETAS)
cbar.set_ticklabels([f"{b:g}" for b in BETAS])

fname = os.path.join(PLOTS, "uncoupled_size_comparison.pdf")
plt.savefig(fname)
print(f"Saved {fname}")
plt.close()

print("Done.")
