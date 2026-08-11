#!/usr/bin/env python3
"""Check: analytically continue the Paper ISO correlations at several beta and
plot the real-time correlation C(t).  One curve per beta; separate Re and Im panels.

Run from the repository root:
    python3 spinDMFT/check_analytic_continuation_betas.py
"""
import os
import sys

import h5py
import numpy as np

# import the continuation machinery from the repo-root script
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analytic_continuation import (
    load_correlations, sanitise_errors, maxent_spectrum,
    real_time_correlation, frequency_grid, _find_matplotlibrc,
)

BETAS = [0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 30.0]
DATADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Data", "Paper")
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Plots")

WMAX, NW = 25.0, 800
TMAX, NT = 10.0, 500          # common real-time axis for all beta

t = np.linspace(0.0, TMAX, NT)
# sinh mesh dense near omega=0 -- resolves the low-T quasi-elastic peak
w = frequency_grid(WMAX, NW, mesh="sinh", conc=4.0)

results = []
for beta in BETAS:
    tag = f"{beta:g}"
    path = os.path.join(DATADIR, f"spinmodel=ISO__beta={tag}.hdf5")
    b, tau, comps = load_correlations(path, part="re")
    comp = comps[0]                       # single local correlation
    g = 0.5 * (comp["g"] + comp["g"][::-1])   # enforce tau <-> beta-tau symmetry
    err, _ = sanitise_errors(g, comp["std"])
    A, g_fit = maxent_spectrum(b, tau, g, err, w, alpha="chi2kink")
    C = real_time_correlation(b, w, A, t)
    chi = np.sqrt(np.mean(((g_fit - g) / err) ** 2))
    results.append((beta, C))
    print(f"beta={beta:>4}:  C(0)={C[0].real:.5f} (g(0)={g[0]:.5f}), "
          f"Im C(0)={C[0].imag:+.1e},  fit chi/pt={chi:.3f}")

# ------------------------------------------------------------------ plotting
import matplotlib
matplotlib.use("Agg")
mplrc = _find_matplotlibrc(os.path.dirname(os.path.abspath(__file__)))
if mplrc:
    matplotlib.rc_file(mplrc)
import matplotlib.pyplot as plt

cmap = plt.get_cmap("viridis")
colors = [cmap(i / (len(BETAS) - 1)) for i in range(len(BETAS))]

fig, ax = plt.subplots(1, 2, figsize=(11, 4.0))
for (beta, C), col in zip(results, colors):
    ax[0].plot(t, C.real, color=col, label=rf"$\beta={beta:g}$")
    ax[1].plot(t, C.imag, color=col, label=rf"$\beta={beta:g}$")

ax[0].set_xlabel(r"$t$")
ax[0].set_ylabel(r"$\mathrm{Re}\,C(t)$")
ax[1].set_xlabel(r"$t$")
ax[1].set_ylabel(r"$\mathrm{Im}\,C(t)$")
for a in ax:
    a.set_xlim(t[0], t[-1])
ax[0].legend()
fig.tight_layout()

os.makedirs(OUTDIR, exist_ok=True)
out = os.path.join(OUTDIR, "analytic_continuation_ISO_vs_beta.pdf")
fig.savefig(out)
print("plot ->", os.path.relpath(out, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
