import glob
import os
import re

import h5py as h5
import matplotlib.pyplot as plt
import numpy as np

DATA_DIR = "Data/Cheb_random"
BETA_ARRAY = [2.5]

_base = plt.rcParams["font.size"]
plt.rcParams.update({
    "font.size": _base + 2,
    "axes.titlesize": _base * 1.2 + 2,
    "axes.labelsize": _base * 1.0 + 4,
    "xtick.labelsize": plt.rcParams["xtick.labelsize"] + 3,
    "ytick.labelsize": plt.rcParams["ytick.labelsize"] + 3,
    "legend.fontsize": _base * 0.833 + 2,
})


def find_files():
    files = {}
    for path in glob.glob(f"{DATA_DIR}/ISO__Random__N=*__beta=*__numConfigs=*.hdf5"):
        m = re.search(r"N=(\d+)__beta=([\d.]+)", path)
        N, beta = int(m.group(1)), float(m.group(2))
        files[(N, beta)] = path
    return files


def load(path):
    f = h5.File(path, "r")
    params = f["parameters"].attrs
    t = np.linspace(0.0, params["Tmax"], params["num_TimePoints"])
    re = f["results"]["Re_correlation"][0]
    re_err = f["results"]["Re_stddev"][0]
    im = f["results"]["Im_correlation"][0]
    im_err = f["results"]["Im_stddev"][0]
    return t, re, re_err, im, im_err


def load_spindmft(beta):
    path = f"{DATA_DIR}/spinmodel=ISO__beta={beta}__forcomparison.hdf5"
    if not os.path.exists(path):
        return None
    with h5.File(path, "r") as f:
        params = f["parameters"].attrs
        t = np.linspace(0.0, params["Tmax"], params["num_RealTimePoints"])
        # The tau=0 slice X^{xx}(t, 0) is the real-time correlator.
        c_real = f["results"]["Re_correlation"][:, 0, 0]
        c_imag = f["results"]["Im_correlation"][:, 0, 0]
        c_real_err = f["runtimedata"]["Re_correlation_sample_stds"][:, 0, 0]
        c_imag_err = f["runtimedata"]["Im_correlation_sample_stds"][:, 0, 0]
    return t, c_real, c_real_err, c_imag, c_imag_err


def extrapolate_1_over_N(Ns, values, errors):
    """Weighted least-squares fit of y = C + a/N at each timepoint, returns C(t), C_err(t)."""
    x = np.array([1.0 / N for N in Ns])
    values = np.array(values)
    errors = np.where(np.array(errors) == 0, 1e-12, np.array(errors))
    w = 1.0 / errors**2
    S = np.sum(w, axis=0)
    Sx = np.sum(w * x[:, None], axis=0)
    Sxx = np.sum(w * (x**2)[:, None], axis=0)
    Sy = np.sum(w * values, axis=0)
    Sxy = np.sum(w * x[:, None] * values, axis=0)
    Delta = S * Sxx - Sx**2
    C = (Sxx * Sy - Sx * Sxy) / Delta
    C_err = np.sqrt(Sxx / Delta)
    return C, C_err


files = find_files()
N_array = sorted({N for N, _ in files})

color_cycle = [c for c in plt.rcParams["axes.prop_cycle"].by_key()["color"] if c != "red"]
colors = {N: color_cycle[i % len(color_cycle)] for i, N in enumerate(N_array)}

PARTS = [
    ("re", 1, 2, r"Re $g^{xx}(t)$"),
    ("im", 3, 4, r"Im $g^{xx}(t)$"),
]

fig, axes = plt.subplots(
    len(PARTS), len(BETA_ARRAY), figsize=(5 * len(BETA_ARRAY), 4.5 * len(PARTS)),
    sharex="col", sharey="row", squeeze=False,
)

for row, (kind, value_idx, err_idx, ylabel) in enumerate(PARTS):
    for col, beta in enumerate(BETA_ARRAY):
        ax = axes[row, col]
        t_common = None
        fit_values, fit_errors, fit_Ns = [], [], []
        for N in N_array:
            key = (N, beta)
            if key not in files:
                continue
            data = load(files[key])
            t_common = data[0]
            ax.errorbar(
                data[0], data[value_idx], yerr=data[err_idx], errorevery=5, capsize=2,
                color=colors[N], linestyle="-", label=rf"$N={N}$",
            )
            fit_values.append(data[value_idx])
            fit_errors.append(data[err_idx])
            fit_Ns.append(N)

        if len(fit_Ns) >= 3:
            C, C_err = extrapolate_1_over_N(fit_Ns, fit_values, fit_errors)
            ax.errorbar(
                t_common, C, yerr=C_err, errorevery=5, capsize=2,
                color="black", linestyle="--", zorder=11, label=r"$N=\infty$",
            )

        spindmft = load_spindmft(beta)
        if spindmft is not None:
            t_dmft, c_real, c_real_err, c_imag, c_imag_err = spindmft
            c, c_err = (c_real, c_real_err) if kind == "re" else (c_imag, c_imag_err)
            ax.errorbar(
                t_dmft, c, yerr=c_err, errorevery=10, capsize=2,
                color="red", linestyle="-", linewidth=2.5, zorder=0,
                label="spinDMFT",
            )

        ax.set_xlim(0, 10)
        if row == len(PARTS) - 1:
            ax.set_xlabel(r"$t J_Q$")
        if col == 0:
            ax.set_ylabel(ylabel)

axes[0, 0].legend()

fig.tight_layout()
fig.subplots_adjust(hspace=0)
fig.savefig("Plots/Plot_realtime_Cheb_random_beta=2.5.pdf")
fig.savefig("Plots/Plot_realtime_Cheb_random_beta=2.5.png", dpi=300)
