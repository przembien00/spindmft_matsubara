#!/usr/bin/env python3
r"""Numerical analytic continuation of imaginary-time spin correlations to real time.

Both the ``spinDMFT`` and ``CspinDMFT_finite_T`` algorithms store the finite-temperature
imaginary-time spin correlation

    g^{ab}_{ij}(tau) = <S^a_i(tau) S^b_j(0)>,     tau in [0, beta]

on an HDF5 file.  This script performs the (ill-posed) analytic continuation

    g(tau) = int_0^inf dw  K(tau, w) A(w)                      (imaginary time)
    C(t)   = int_0^inf dw  K(t/beta -> i(t/beta), beta w) A(w) (real time)

using Maximum-Entropy (MaxEnt) as implemented in the pre-made library ``ana_cont``
(J. Kaufmann & K. Held, Comput. Phys. Commun. 282, 108519 (2023)).

Pipeline
--------
1.  Read g(tau) (real part by default) and its Monte-Carlo standard deviations.
2.  MaxEnt with the *bosonic* imaginary-time kernel  ->  non-negative spectral
    function  A(w)  on the real-frequency axis.  The alpha hyper-parameter is
    chosen automatically with the chi2-kink method.
3.  Reconstruct the real-time correlation  C(t)  as the analytic continuation
    tau -> i t of the *same* kernel the library fits, so that C(0) = g(0) exactly
    and Im C(0) = 0 by construction.

The bosonic kernel assumes a correlation that is real, symmetric g(tau)=g(beta-tau)
and non-negative (A(w) >= 0).  This is exactly the case for the on-site
autocorrelation g^{aa}_{ii}; off-site / off-diagonal components that change sign
are continued too but a warning is printed since the positive-spectrum assumption
may not hold there.

Real-frequency / real-time conventions
--------------------------------------
The library's dimensionless bosonic imaginary-time kernel (x = beta*w, s = tau/beta) is

    K(s, x) = 0.5 * x * [ e^{-x s} + e^{-x (1-s)} ] / (1 - e^{-x}) ,   K(0-limit)=1 .

Continuing s -> i u  with  u = t/beta  gives the real-time kernel

    K_R(u, x) = 0.5 * x * [ e^{-i x u} + e^{-x} e^{+i x u} ] / (1 - e^{-x})

and the real-time correlation

    C(t) = (1/beta) * int_0^inf dw  K_R(t/beta, beta w)  A(w) .

Re C(t) is the symmetrized (Kubo) correlation; Im C(t) the response part.

Usage
-----
    python3 analytic_continuation.py <file.hdf5> [options]

    # continue every stored component, all default settings
    python3 analytic_continuation.py CspinDMFT_finite_T/Data/.../beta=2.5.hdf5

    # a single component with a custom frequency window
    python3 analytic_continuation.py file.hdf5 --component 0-0 --wmax 30 --tmax 15

Options: run  python3 analytic_continuation.py --help .

Output
------
* An HDF5 file  <input>__realtime.hdf5  (next to the input) holding, per component,
  the real-frequency axis + A(w) and the real-time axis + Re/Im C(t).
* A plot per component in a  Plots/  folder beside the input file (unless --noplot).
"""

import argparse
import os
import sys
import warnings

import h5py
import numpy as np
from scipy.integrate import trapezoid

try:
    from ana_cont import continuation as cont
except ImportError:
    sys.exit(
        "ana_cont is required.  Install it with:\n"
        "    pip3 install --user Cython\n"
        "    pip3 install --user ana_cont\n"
    )


# --------------------------------------------------------------------------- #
#  I/O: read the imaginary-time correlations from either algorithm's HDF5
# --------------------------------------------------------------------------- #
def load_correlations(path, part="re"):
    """Return (beta, tau, components) for a spinDMFT or CspinDMFT_finite_T file.

    ``components`` is a list of dicts with keys:
        label : str          human-readable id, e.g. '0-0' or '0-0.c1'
        g     : (ntau,) array imaginary-time correlation g(tau)
        std   : (ntau,) array Monte-Carlo standard deviation (>=0)
    ``part`` selects the real ('re') or imaginary ('im') part of g(tau).
    """
    prefix = {"re": "Re", "im": "Im"}[part]
    with h5py.File(path, "r") as h:
        beta = float(h["parameters"].attrs["beta"])
        res = h["results"]
        corr_name = f"{prefix}_correlation"
        node = res[corr_name]

        components = []
        if isinstance(node, h5py.Dataset):
            # ---- spinDMFT: single local correlation, shape (1, ntau) ----
            g_all = np.asarray(node)
            # newer (Paper) files split stds by part; older files use one dataset
            std_node = h.get(f"runtimedata/{prefix}_correlation_sample_stds")
            if std_node is None:
                std_node = h.get("runtimedata/correlation_sample_stds")
            std_all = np.asarray(std_node) if std_node is not None else np.zeros_like(g_all)
            for c in range(g_all.shape[0]):
                components.append(dict(label="local" if g_all.shape[0] == 1 else f"c{c}",
                                       g=g_all[c].astype(float),
                                       std=std_all[c].astype(float)))
        else:
            # ---- CspinDMFT_finite_T: group of site-pair datasets 'i-j' ----
            std_grp = h.get(f"runtimedata/{prefix}_correlation_sample_stds")
            for key in sorted(node.keys(), key=_pair_key):
                g_all = np.asarray(node[key])              # (ncomp, ntau)
                if std_grp is not None and key in std_grp:
                    std_all = np.asarray(std_grp[key])
                else:
                    std_all = np.zeros_like(g_all)
                ncomp = g_all.shape[0]
                for c in range(ncomp):
                    label = key if ncomp == 1 else f"{key}.c{c}"
                    components.append(dict(label=label,
                                           g=g_all[c].astype(float),
                                           std=std_all[c].astype(float)))

    ntau = components[0]["g"].size
    tau = np.linspace(0.0, beta, ntau)
    return beta, tau, components


def _pair_key(s):
    """Sort key so '0-1' < '0-10' < '2-3' numerically."""
    try:
        i, j = s.split("-")
        return (int(i), int(j))
    except ValueError:
        return (0, s)


# --------------------------------------------------------------------------- #
#  Real-frequency mesh
# --------------------------------------------------------------------------- #
def frequency_grid(wmax, nw, mesh="sinh", conc=4.0):
    """Real-frequency grid on [0, wmax].

    Spin correlations here have all their spectral weight near omega = 0 (the
    quasi-elastic / slow mode), and at low temperature that peak is only ~1/beta
    wide.  A uniform grid wastes resolution: the default 'sinh' mesh places points
    densely near 0 and sparsely at large omega while still covering [0, wmax].

    ``conc`` controls the concentration (larger -> more points near 0).  The mesh
    must start exactly at 0 (the bosonic kernel's de l'Hospital point).
    """
    if mesh == "linear":
        return np.linspace(0.0, wmax, nw)
    if mesh == "sinh":
        u = np.linspace(0.0, 1.0, nw)
        return wmax * np.sinh(conc * u) / np.sinh(conc)
    raise ValueError(f"unknown mesh '{mesh}' (use 'sinh' or 'linear')")


# --------------------------------------------------------------------------- #
#  Analytic continuation
# --------------------------------------------------------------------------- #
def sanitise_errors(g, std, err_floor=1e-4, default_rel=1e-3):
    """Return a usable per-tau error array.

    MaxEnt needs the error bars to reflect the real precision of g(tau).  Two
    pathologies appear in the raw output and are corrected here:

    * exact zeros at some tau (would divide by zero)  -> floored;
    * the spinDMFT ``correlation_sample_stds`` dataset is not a normalised
      standard deviation (values ~1e5 while g ~ 1); when std is absurdly large
      compared with g it is discarded in favour of a flat ``default_rel`` error.
    """
    scale = max(abs(g[0]), 1e-12)
    std = np.asarray(std, dtype=float)
    if not np.any(std > 0) or np.median(std) > np.ptp(g) + scale:
        return np.full_like(g, default_rel * scale), True
    return np.clip(std, err_floor * scale, None), False


def maxent_spectrum(beta, tau, g, err, w, alpha="chi2kink"):
    """Run bosonic MaxEnt; return the spectral function A(w) and the reconstructed g.

    If the requested alpha-determination fails to find a kink (e.g. very clean or
    over-fitted data), fall back to Bryan's alpha-averaged solution.
    """
    # Flat default model.  The bosonic ("positive-negative") entropy lets the norm
    # float, so the exact model normalisation is irrelevant; keep it unit-normalised.
    model = np.ones_like(w)
    model /= trapezoid(model, w)

    problem = cont.AnalyticContinuationProblem(
        im_axis=tau, re_axis=w, im_data=g, kernel_mode="time_bosonic", beta=beta
    )

    def _solve(method):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol, _ = problem.solve(
                method="maxent_svd", model=model, stdev=err,
                alpha_determination=method, optimizer="newton",
            )
        return sol

    try:
        sol = _solve(alpha)
    except RuntimeError as exc:
        print(f"           alpha='{alpha}' failed ({exc}); falling back to 'bryan'.")
        sol = _solve("bryan")
    return sol.A_opt, sol.backtransform


def real_time_correlation(beta, w, A, t):
    """C(t) via analytic continuation tau -> i t of the bosonic kernel.

    C(t) = (1/beta) int dw  0.5 x [e^{-i x u} + e^{-x} e^{+i x u}]/(1-e^{-x}) A(w),
    with x = beta w and u = t/beta.  Returns a complex array of shape t.shape.
    """
    x = beta * w
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = 1.0 - np.exp(-x)
    C = np.empty(t.shape, dtype=complex)
    for k, tt in enumerate(t):
        u = tt / beta
        K = 0.5 * x * (np.exp(-1j * x * u) + np.exp(-x) * np.exp(1j * x * u)) / denom
        K[x == 0.0] = 1.0  # de l'Hospital limit, matches the library
        C[k] = trapezoid(K * A, w) / beta
    return C


# --------------------------------------------------------------------------- #
#  Plotting (follows the repository conventions)
# --------------------------------------------------------------------------- #
def _find_matplotlibrc(start):
    """Walk up from ``start`` to locate the algorithm's matplotlibrc, if any."""
    d = os.path.abspath(start)
    while True:
        cand = os.path.join(d, "matplotlibrc")
        if os.path.isfile(cand):
            return cand
        parent = os.path.dirname(d)
        if parent == d:
            return None
        d = parent


def plot_component(label, beta, tau, g, g_fit, w, A, t, C, outdir, subscript, mplrc=None):
    import matplotlib
    matplotlib.use("Agg")
    if mplrc:
        matplotlib.rc_file(mplrc)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))

    # (1) imaginary-time data + MaxEnt reconstruction, tau normalised to [0,1]
    ax[0].plot(tau / beta, g, "o", ms=3, label="data")
    ax[0].plot(tau / beta, g_fit, "-", lw=1.5, label="MaxEnt fit")
    ax[0].set_xlabel(r"$\tau/\beta$")
    ax[0].set_ylabel(rf"$g^{{{subscript}}}(\tau)$")
    ax[0].set_xlim(0.0, 1.0)
    ax[0].legend()

    # (2) spectral function A(w)
    ax[1].plot(w, A, "-", lw=1.5)
    ax[1].set_xlabel(r"$\omega$")
    ax[1].set_ylabel(r"$A(\omega)$")
    ax[1].set_xlim(w[0], w[-1])

    # (3) real-time correlation
    ax[2].plot(t, C.real, "-", lw=1.5, label=r"$\mathrm{Re}\,C(t)$")
    ax[2].plot(t, C.imag, "--", lw=1.5, label=r"$\mathrm{Im}\,C(t)$")
    ax[2].set_xlabel(r"$t$")
    ax[2].set_ylabel(rf"$C^{{{subscript}}}(t)$")
    ax[2].set_xlim(t[0], t[-1])
    ax[2].legend()

    fig.tight_layout()
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"analytic_continuation_{label}.pdf")
    fig.savefig(out)
    plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
#  Main
# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(
        description="Analytic continuation of imaginary-time spin correlations to real time.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("file", help="input HDF5 file (spinDMFT or CspinDMFT_finite_T)")
    p.add_argument("--component", default="all",
                   help="component label to continue (e.g. '0-0'), or 'all'")
    p.add_argument("--part", choices=["re", "im"], default="re",
                   help="real or imaginary part of g(tau) to continue")
    p.add_argument("--wmax", type=float, default=25.0, help="max real frequency omega")
    p.add_argument("--nw", type=int, default=800, help="number of frequency points")
    p.add_argument("--mesh", choices=["sinh", "linear"], default="sinh",
                   help="frequency mesh: 'sinh' is dense near omega=0 (better for "
                        "the omega=0-peaked spin spectra, especially at low T)")
    p.add_argument("--mesh-conc", type=float, default=4.0,
                   help="sinh-mesh concentration (larger -> denser near omega=0)")
    p.add_argument("--tmax", type=float, default=None,
                   help="max real time (default: 5*beta)")
    p.add_argument("--nt", type=int, default=400, help="number of real-time points")
    p.add_argument("--alpha", default="chi2kink",
                   choices=["chi2kink", "classic", "bryan", "historic"],
                   help="MaxEnt alpha-determination method")
    p.add_argument("--outdir", default=None,
                   help="directory for plots (default: <input dir>/Plots)")
    p.add_argument("--noplot", action="store_true", help="skip plotting")
    args = p.parse_args()

    if not os.path.isfile(args.file):
        sys.exit(f"file not found: {args.file}")

    beta, tau, components = load_correlations(args.file, part=args.part)
    if args.component != "all":
        components = [c for c in components if c["label"] == args.component]
        if not components:
            sys.exit(f"component '{args.component}' not found in file")

    w = frequency_grid(args.wmax, args.nw, mesh=args.mesh, conc=args.mesh_conc)
    tmax = args.tmax if args.tmax is not None else 5.0 * beta
    t = np.linspace(0.0, tmax, args.nt)

    base = os.path.splitext(os.path.basename(args.file))[0]
    indir = os.path.dirname(os.path.abspath(args.file))
    outdir = args.outdir or os.path.join(indir, "Plots")
    out_h5 = os.path.join(indir, base + "__realtime.hdf5")
    mplrc = None if args.noplot else _find_matplotlibrc(indir)

    print(f"beta = {beta},  ntau = {tau.size},  components = {len(components)}")
    with h5py.File(out_h5, "w") as hout:
        hout.attrs["source_file"] = os.path.abspath(args.file)
        hout.attrs["beta"] = beta
        hout.attrs["part"] = args.part
        hout.attrs["method"] = "ana_cont MaxEnt (time_bosonic kernel)"
        hout.create_dataset("real_frequency", data=w)
        hout.create_dataset("real_time", data=t)

        for comp in components:
            label = comp["label"]
            g, std = comp["g"], comp["std"]

            # sanity checks for the bosonic (positive, symmetric) assumption
            asym = np.max(np.abs(g - g[::-1])) / max(np.max(np.abs(g)), 1e-12)
            if asym > 1e-3:
                print(f"  [{label}] WARNING: g(tau) not symmetric about beta/2 "
                      f"(rel. asym {asym:.1e}); symmetrising.")
                g = 0.5 * (g + g[::-1])
            if g.min() < -1e-6 * abs(g[0]):
                print(f"  [{label}] WARNING: g(tau) changes sign; the bosonic "
                      f"positive-spectrum assumption may not hold.")

            err, replaced = sanitise_errors(g, std)
            if replaced:
                print(f"  [{label}] NOTE: sample_stds unusable; using flat 1e-3 "
                      f"relative error for MaxEnt.")
            A, g_fit = maxent_spectrum(beta, tau, g, err, w, alpha=args.alpha)
            C = real_time_correlation(beta, w, A, t)

            chi = np.sqrt(np.mean(((g_fit - g) / err) ** 2))
            print(f"  [{label}]  C(0)={C[0].real:+.5f} (g(0)={g[0]:+.5f}),  "
                  f"Im C(0)={C[0].imag:+.1e},  fit chi/pt={chi:.2e}")

            grp = hout.create_group(label)
            grp.create_dataset("spectral_function", data=A)
            grp.create_dataset("C_real", data=C.real)
            grp.create_dataset("C_imag", data=C.imag)
            grp.create_dataset("g_tau", data=g)
            grp.create_dataset("g_tau_fit", data=g_fit)

            if not args.noplot:
                subscript = "aa_{" + label + "}"
                path = plot_component(label, beta, tau, g, g_fit, w, A, t, C,
                                      outdir, subscript, mplrc=mplrc)
                print(f"           plot -> {os.path.relpath(path)}")

    print(f"results -> {os.path.relpath(out_h5)}")


if __name__ == "__main__":
    main()
