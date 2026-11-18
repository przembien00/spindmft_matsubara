"""Compare cubic-lattice real-time correlations and stored uncertainties."""

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
BETAS = (0.2, 1.5, 2.5)
TMAX = 10
CASES = {
    "AFM_B0": ("AFM/zero", "Jplus", "0"),
    "AFM_B0p5": ("AFM/hz0p5", "Jplus", "0.5"),
    "FM_B0": ("FM/zero", "Jminus", "0"),
    "FM_B0p5": ("FM/hz0p5", "Jminus", "0.5"),
}
# Representatives from Cubic3D_NSpin.py, indexed as in the Keldysh output.
PAIRS = {
    "onsite": {2: "1-1", 4: "1-1", 7: "1-1"},
    "NN": {2: "1-2", 4: "1-2", 7: "1-2"},
    "NNN": {4: "1-4", 7: "2-4"},
}
QMC_PAIRS = {"onsite": "0-0", "NN": "1-0", "NNN": "33-0"}
COLORS = {2: "#004dff", 4: "#ff7800", 7: "#00a824"}
PARTS = ("Re", "Im")


def select_runs(directory):
    runs = {beta: {} for beta in BETAS}
    for path in sorted(directory.glob("*.hdf5")):
        with h5py.File(path, "r") as source:
            params = source["parameters"].attrs
            beta = float(params["beta"])
            if beta not in runs or "results" not in source:
                continue
            size = int(params["num_Spins"])
            samples = int(params["num_Samples"])
            if size not in runs[beta] or samples > runs[beta][size][0]:
                runs[beta][size] = (samples, path)
    return runs


def select_single_site_runs(field):
    runs = {beta: {} for beta in BETAS}
    for path in sorted((ROOT / "Data/Single-site").glob("*.hdf5")):
        with h5py.File(path, "r") as source:
            params = source["parameters"].attrs
            beta = float(params["beta"])
            jl = float(params["JL"])
            matches_coupling = jl == 0 if field == "0" else jl < 0
            if beta in runs and matches_coupling and float(params.get("h_abs", 0)) == float(field):
                runs[beta][jl] = path
    return runs


def load_keldysh(path, pair, component):
    with h5py.File(path, "r") as source:
        params = source["parameters"].attrs
        results = source["results"]
        index = list(results["correlation_direction_labels"].asstr()[:]).index(component)
        time = results.attrs.get("real_times")
        if time is None:
            time = np.linspace(0, params["Tmax"], params["num_RealTimePoints"])
        pair_suffix = f"/{pair}" if pair else ""
        keep = time <= TMAX
        series = []
        for part in PARTS:
            values = results[f"{part}_correlation{pair_suffix}"][keep, index, 0]
            errors = source[f"runtimedata/{part}_correlation_sample_stds{pair_suffix}"][keep, index, 0]
            assert np.isfinite(values).all() and np.isfinite(errors).all()
            assert (errors >= 0).all()
            series.append((values, values - errors, values + errors))
    return time[keep], series


def load_qmc(group, pair, component):
    key = f"Re_g/{pair}"
    if key not in group:
        return None
    index = group[key].attrs["spin_components"].split(",").index(component)
    time = group["t"][:]
    # Retain the point just beyond TMAX so the line reaches the axis boundary.
    stop = np.searchsorted(time, TMAX, side="right") + 1
    series = []
    for part in PARTS:
        values = group[f"{part}_g/{pair}"][index, :stop]
        low = group[f"{part}_g_band_low/{pair}"][index, :stop]
        high = group[f"{part}_g_band_high/{pair}"][index, :stop]
        assert np.isfinite([values, low, high]).all() and (low <= high).all()
        # Conjugate to the Keldysh ordering, reversing the imaginary bounds.
        series.append((values, low, high) if part == "Re" else (-values, -high, -low))
    return time[:stop], series


def draw_series(axes, data, color, *, qmc=False):
    time, series = data
    for ax, (values, low, high) in zip(axes, series):
        ax.fill_between(time, low, high, color=color, alpha=0.18, linewidth=0,
                        zorder=1 if qmc else 2)
        (handle,) = ax.plot(time, values, color=color, linestyle="-" if qmc else "--",
                            linewidth=1.1 if qmc else 1.35, zorder=3 if qmc else 4)
    return handle


def plot_case(name, runs, single_site, source, shell, component, output):
    fig, axes = plt.subplots(2, 3, figsize=(12, 6), sharex=True, sharey="row")
    handles = {}
    for column, beta in enumerate(BETAS):
        panels = axes[:, column]
        sizes = sorted(runs[beta].keys() & PAIRS[shell].keys())
        for size in sizes:
            _, path = runs[beta][size]
            data = load_keldysh(path, PAIRS[shell][size], component)
            handles[f"CspinDMFT $N={size}$"] = draw_series(panels, data, COLORS[size])
        for jl, path in sorted(single_site.get(beta, {}).items()):
            label = rf"Single-site $J_L={jl:g}$"
            color = "#008ba3" if jl else "#a66b00"
            data = load_keldysh(path, None, component)
            handles[label] = draw_series(panels, data, color)
        qmc = load_qmc(source[f"beta={beta:g}"], QMC_PAIRS[shell], component)
        if qmc is not None:
            handles["QMC"] = draw_series(panels, qmc, "black", qmc=True)
        if not sizes and not single_site.get(beta) and qmc is None:
            for ax in panels:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
        panels[0].text(0.96, 0.93, rf"$\beta={beta:g}$",
                       transform=panels[0].transAxes, ha="right", va="top")
        panels[1].set_xlabel(r"$t$")
        print(name, shell, component, beta, sizes, sorted(single_site.get(beta, {})),
              "QMC" if qmc is not None else "")
    axes[0, 0].set_xlim(0, TMAX)
    pair_label = "ii" if shell == "onsite" else "ij"
    shell_label = "" if shell == "onsite" else f" ({shell})"
    for ax, part in zip(axes[:, 0], PARTS):
        ax.set_ylabel(rf"$\mathrm{{{part}}}\, g^{{{component}}}_{{{pair_label}}}(t)$" + shell_label)
    labels = sorted(handles, key=lambda label: (label != "QMC", label))
    if "QMC" not in handles:
        fig.text(0.5, 0.065, "QMC data unavailable for this neighbor shell",
                 ha="center", fontsize="small")
    fig.legend([handles[label] for label in labels], labels,
               loc="lower center", ncol=len(labels), frameon=False)
    fig.tight_layout(rect=(0, 0.095, 1, 1))
    stem = name if (shell, component) == ("onsite", "xx") else f"{name}_{shell}_{component}"
    fig.savefig(output / f"{stem}_qmc_comparison.pdf")
    plt.close(fig)


def main():
    output = ROOT / "Plots/Cubic_QMC_comparison"
    output.mkdir(parents=True, exist_ok=True)
    plt.style.use(ROOT.parent / "spinDMFT_Keldysh/matplotlibrc")
    plt.rcParams["figure.autolayout"] = False
    case_folders = {"AFM_B0": "AFM", "AFM_B0p5": "AFM_B",
                    "FM_B0": "FM", "FM_B0p5": "FM_B"}
    for name, (directory, coupling, field) in CASES.items():
        runs = select_runs(ROOT / "Data/Cubic_beta_sweep" / directory)
        qmc_path = ROOT / "Data/QMC" / f"realtime_{coupling}_hz{field}_N32768_all_betas.hdf5"
        components = ("xx", "xy", "zz") if field != "0" else ("xx",)
        single_site = {} if name == "AFM_B0p5" else select_single_site_runs(field)
        with h5py.File(qmc_path, "r") as source:
            for shell in PAIRS:
                shell_folder = "local" if shell == "onsite" else shell
                destination = output / case_folders[name] / shell_folder
                destination.mkdir(parents=True, exist_ok=True)
                for component in components:
                    local_runs = single_site if shell == "onsite" else {}
                    plot_case(name, runs, local_runs, source, shell, component, destination)


if __name__ == "__main__":
    main()
