#!/usr/bin/env python3
"""Compare spinDMFT (N=2,4,9) vs Chebyshev (N=20) for AFM and FM.

Six figures: {AFM,FM} × {local, NN, NNN}.
Color = beta J (plasma_r: yellow=small, purple=large).
Line style = N; solid = Chebyshev, dashed/dashdot/dotted = CspinDMFT.
"""

import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
import matplotlib.colorbar as mcolorbar

ROOT  = os.path.dirname(os.path.abspath(__file__))
DATA  = os.path.join(ROOT, "Data")
PLOTS = os.path.join(ROOT, "Plots", "size_comparison")
os.makedirs(PLOTS, exist_ok=True)

BETAS = [0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0][::2]

# solid = Chebyshev (N=24), non-solid = CspinDMFT
N_LSTYLE = {2: "--", 4: "-.", 9: ":", 24: "-"}

CMAP      = plt.get_cmap("plasma_r")
BETA_NORM = mcolors.Normalize(vmin=min(BETAS), vmax=max(BETAS))

def beta_color(beta: float):
    return CMAP(BETA_NORM(beta))

# ── DMFT loaders ──────────────────────────────────────────────────────────────
def dmft_path(mag_type: str, N: int, beta: float) -> str:
    J  = 0.5  if mag_type == "AFM" else -0.5
    bstr = f"{beta:g}"
    return os.path.join(DATA, mag_type,
        f"spinmodel=ISO__config=Square_2D_N={N}_NN_J={J}__beta={bstr}.hdf5")

def load_dmft(path: str, key: str):
    with h5py.File(path, "r") as f:
        delta_t = float(f["parameters"].attrs["delta_t"])
        ntp     = f["results/Re_correlation/0-0"].shape[1]
        data    = f["results/Re_correlation"][key][()][0].astype(float)
        stds    = f["runtimedata/Re_correlation_sample_stds"][key][()][0].astype(float)
        n_samp  = float(f["parameters"].attrs["num_Samples"])
    tau = np.arange(ntp) * delta_t / (delta_t * (ntp - 1))  # normalise to [0,1]
    return tau, data, stds / np.sqrt(n_samp)

# ── Chebyshev loaders ─────────────────────────────────────────────────────────
def cheb_path(mag_type: str, beta: float) -> str:
    rescale = 0.5 if mag_type == "AFM" else -0.5
    bstr = f"{beta:g}"
    return os.path.join(DATA, "Chebyshev",
        f"ISO__Square_NN_PBC_N=24__sites=0-1-7__beta={bstr}__rescale={rescale}.hdf5")

def load_cheb(path: str, key: str):
    with h5py.File(path, "r") as f:
        data = f["results/Re_correlation"][key][()][0].astype(float)
        std  = f["results/Re_stds"][key][()][0].astype(float)
    data = np.concatenate([data, data[-2::-1]])
    std  = np.concatenate([std,  std[-2::-1]])
    tau  = np.linspace(0, 1, len(data))
    return tau, data, std

# ── correlation key specifications ────────────────────────────────────────────
# corr_spec[corr_name] = dict(dmft_keys, cheb_key, ylabel, tag)
# dmft_keys: N → key string (None = not available)
CORR_SPECS = {
    "local": dict(
        dmft_keys={2: "0-0", 4: "0-0", 9: "0-0"},
        cheb_key="0-0",
        ylabel=r"$g^{xx}_\mathrm{local}(\tau/\beta)$",
        tag="local",
    ),
    "NN": dict(
        dmft_keys={2: "0-1", 4: "0-1", 9: "0-1"},
        cheb_key="1-0",
        ylabel=r"$g^{xx}_\mathrm{NN}(\tau/\beta)$",
        tag="NN",
    ),
    "NNN": dict(
        dmft_keys={2: None, 4: "0-3", 9: "0-4"},
        cheb_key="7-0",
        ylabel=r"$g^{xx}_\mathrm{NNN}(\tau/\beta)$",
        tag="NNN",
    ),
}

# ── main plotting loop ─────────────────────────────────────────────────────────
for mag_type in ("AFM", "FM"):
    for corr_name, spec in CORR_SPECS.items():
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ymin, ymax = np.inf, -np.inf

        for beta in BETAS:
            # ── DMFT N=2,4,9 ──
            for N in (2, 4, 9):
                dmft_key = spec["dmft_keys"][N]
                if dmft_key is None:
                    continue
                path = dmft_path(mag_type, N, beta)
                if not os.path.exists(path):
                    continue
                try:
                    tau, data, err = load_dmft(path, dmft_key)
                except (KeyError, OSError):
                    continue
                color = beta_color(beta)
                ax.plot(tau, data,
                        color=color, lw=1.2, ls=N_LSTYLE[N])
                ax.fill_between(tau, data - err, data + err,
                                color=color, alpha=0.12)
                ymin = min(ymin, np.min(data))
                ymax = max(ymax, np.max(data))

            # ── Chebyshev N=20 ──
            path = cheb_path(mag_type, beta)
            if not os.path.exists(path):
                continue
            try:
                tau, data, err = load_cheb(path, spec["cheb_key"])
            except (KeyError, OSError):
                continue
            color = beta_color(beta)
            ax.plot(tau, data,
                    color=color, lw=1.2, ls=N_LSTYLE[24])
            ax.fill_between(tau, data - err, data + err,
                            color=color, alpha=0.12)
            ymin = min(ymin, np.min(data))
            ymax = max(ymax, np.max(data))

        ax.axhline(0, color="gray", lw=0.5, ls=":")
        ax.set_xlim(0, 1)
        margin = 0.05 * (ymax - ymin)
        ax.set_ylim(ymin - margin, ymax + margin)
        ax.set_xlabel(r"$\tau/\beta$")
        ax.set_ylabel(spec["ylabel"])

        # Legend: line styles for N values
        n_handles = []
        for N, ls in N_LSTYLE.items():
            lbl = f"$N={N}$" + (" (Chebyshev)" if N == 24 else " (spinDMFT)")
            if corr_name == "NNN" and N == 2:
                continue
            n_handles.append(mlines.Line2D([], [], color="k", ls=ls,
                                           lw=1.5, label=lbl))
        ax.legend(handles=n_handles, fontsize=8, loc="best", framealpha=0.85)

        # Colorbar for beta
        sm = cm.ScalarMappable(cmap=CMAP, norm=BETA_NORM)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.02, aspect=30)
        cbar.set_label(r"$\beta J$", fontsize=10)
        cbar.set_ticks(BETAS)
        cbar.set_ticklabels([f"{b:g}" for b in BETAS])

        plt.tight_layout()
        fname = os.path.join(PLOTS, f"size_comparison__{mag_type}__{corr_name}.pdf")
        plt.savefig(fname, bbox_inches="tight")
        print(f"Saved {fname}")
        plt.close()

print("Done.")
