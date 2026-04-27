import os
import re

import h5py as h5
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.style.use("ggplot")


ROOT_IT = "/Users/przembien/Projects/spindmft_imaginary_time/spinDMFT"
ROOT_MA = "/Users/przembien/Projects/spindmft_matsubara/spinDMFT"
PROJECT_IT = "AlgoCompare_IT_JL"
PROJECT_MA = "AlgoCompare_MA_JL"
OUTPUT = os.path.join(ROOT_MA, "Plots", "AlgoCompare_IT_vs_MA_JLpm2_hz05_xx_xy.pdf")

BETAS = [1.0, 1.5, 2.0]
JLS = [2, -2]

COLORS = {
    1.0: "tab:blue",
    1.5: "tab:orange",
    2.0: "tab:green",
}
MARKERS = {
    1.0: "o",
    1.5: "s",
    2.0: "^",
}
LINESTYLES = {
    "IT": "-",
    "MA": "--",
}


def _file_path(root, project, jl, beta):
    beta_str = f"{beta:.2g}"
    return os.path.join(
        root,
        "Data",
        project,
        f"spinmodel=ISO__JL={jl}__beta={beta_str}__h=z_h_abs=0.5.hdf5",
    )


def _load_components(path):
    with h5.File(path, "r") as all_data:
        params = all_data["parameters"]
        disc = np.linspace(0.0, 1.0, params.attrs["num_TimePoints"])
        re_corr = np.array(all_data["results"]["Re_correlation"])
        im_corr = np.array(all_data["results"]["Im_correlation"])
    return disc, re_corr, im_corr


def _find_first_existing(root, project, jl, beta):
    path = _file_path(root, project, jl, beta)
    if os.path.exists(path):
        return path
    return None


def main():
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharex=True)

    for col, jl in enumerate(JLS):
        for beta in BETAS:
            path_it = _find_first_existing(ROOT_IT, PROJECT_IT, jl, beta)
            path_ma = _find_first_existing(ROOT_MA, PROJECT_MA, jl, beta)
            if path_it is None or path_ma is None:
                missing = []
                if path_it is None:
                    missing.append(f"IT beta={beta:g} JL={jl}")
                if path_ma is None:
                    missing.append(f"MA beta={beta:g} JL={jl}")
                raise FileNotFoundError("Missing input file(s): " + ", ".join(missing))

            disc_it, re_it, im_it = _load_components(path_it)
            disc_ma, re_ma, im_ma = _load_components(path_ma)

            xx_it = np.array(re_it[0], dtype=float)
            xx_ma = np.array(re_ma[0], dtype=float)
            xy_it = np.array(im_it[1], dtype=float)
            xy_ma = np.array(im_ma[1], dtype=float)

            color = COLORS[beta]
            marker = MARKERS[beta]
            label_it = rf"IT, $\beta={beta:g}$"
            label_ma = rf"MA, $\beta={beta:g}$"

            axes[0, col].plot(
                disc_it,
                xx_it,
                color=color,
                linestyle=LINESTYLES["IT"],
                marker=marker,
                markevery=max(1, len(disc_it) // 12),
                label=label_it,
            )
            axes[0, col].plot(
                disc_ma,
                xx_ma,
                color=color,
                linestyle=LINESTYLES["MA"],
                marker=marker,
                markevery=max(1, len(disc_ma) // 12),
                label=label_ma,
            )
            axes[1, col].plot(
                disc_it,
                xy_it,
                color=color,
                linestyle=LINESTYLES["IT"],
                marker=marker,
                markevery=max(1, len(disc_it) // 12),
                label=label_it,
            )
            axes[1, col].plot(
                disc_ma,
                xy_ma,
                color=color,
                linestyle=LINESTYLES["MA"],
                marker=marker,
                markevery=max(1, len(disc_ma) // 12),
                label=label_ma,
            )

        axes[0, col].text(
            0.03,
            0.92,
            rf"$J_L={jl}$",
            transform=axes[0, col].transAxes,
        )

    axes[0, 0].set_ylabel(r"$g^{xx}(\tau)$")
    axes[1, 0].set_ylabel(r"$\mathrm{Im}\,g^{xy}(\tau)$")
    axes[1, 0].set_xlabel(r"$\tau/\beta$")
    axes[1, 1].set_xlabel(r"$\tau/\beta$")

    for ax in axes.flat:
        ax.set_xlim(0.0, 1.0)
        ax.margins(x=0.0)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=3,
        fontsize=8,
        frameon=True,
        bbox_to_anchor=(0.5, 1.02),
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))

    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
    fig.savefig(OUTPUT, dpi=1000)


if __name__ == "__main__":
    main()
