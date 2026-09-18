#!/usr/bin/env python3
"""Minimal real-time g^zz comparison for beta=0.5 and beta=2.5."""

from pathlib import Path

import h5py
import matplotlib

# This script only writes files; using a non-GUI backend keeps it reliable on
# headless compute nodes and in batch rendering environments.
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.axisartist import Axes


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data" / "Keldysh_stype_A"
OUTPUT = ROOT / "Plots" / "Keldysh_gzz_temperature_comparison"


def load_correlation(beta: float):
    path = DATA / f"spinmodel=ISO__beta={beta:g}_samples_per_core=1000000.hdf5"
    with h5py.File(path, "r") as h5file:
        labels = h5file["results/correlation_direction_labels"].asstr()[...].tolist()
        # The zero-field ISO state is rotationally invariant; these files store
        # its sole independent diagonal channel as xx.
        component = labels.index("zz") if "zz" in labels else labels.index("xx")
        times = h5file["results"].attrs["real_times"]
        values = h5file["results/Re_correlation"][:, component, 0]
    return times, values


def main() -> None:
    plt.style.use(ROOT / "matplotlibrc")
    plt.rcParams.update({"axes.grid": False, "xtick.minor.visible": False, "ytick.minor.visible": False})

    fig = plt.figure(figsize=(6.4, 4.1))
    ax = fig.add_subplot(axes_class=Axes)
    curves = (
        (0.5, "red", r"$\mathsf{high}\ T$", 2.5, 22, 31),
        (2.5, "blue", r"$\mathsf{low}\ T$", 6.3, 0, 10),
    )
    for beta, color, label, label_time, label_offset_x, label_offset_y in curves:
        times, values = load_correlation(beta)
        ax.plot(times, values, color=color, linewidth=1.8)
        label_index = int(round(label_time / (times[1] - times[0])))
        ax.annotate(
            label,
            xy=(times[label_index], values[label_index]),
            xytext=(label_offset_x, label_offset_y),
            textcoords="offset points",
            color=color,
            ha="center",
            va="bottom",
        )

    ax.set_xlim(0.0, 10.0)
    ax.set_ylim(bottom=0.0)
    ax.axis["top"].set_visible(False)
    ax.axis["right"].set_visible(False)
    ax.axis["bottom"].set_axisline_style("-|>", size=2.5)
    ax.axis["left"].set_axisline_style("-|>", size=2.5)
    ax.axis["bottom"].label.set_visible(False)
    ax.axis["left"].label.set_visible(False)
    label_size = plt.rcParams["axes.labelsize"]
    ax.text(0.99, -0.037, r"$t$", transform=ax.transAxes, ha="right", va="top", fontsize=label_size)
    ax.text(-0.035, 0.99, r"$g^{zz}(t)$", transform=ax.transAxes, ha="right", va="top", fontsize=label_size)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(bottom=False, left=False)
    fig.tight_layout()

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        fig.savefig(OUTPUT.with_suffix(f".{extension}"), dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
