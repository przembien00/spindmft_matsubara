"""Freezing-transition scan of the N=2 uncoupled cluster (J=0.5):
order parameter r = g^{zz}_{01}(beta/2) / g^{zz}_{00}(beta/2) vs beta,
and the freezing-condition check lambda * f(beta_c) = 1 across cluster sizes."""
import glob
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

files = sorted(
    glob.glob("Data/Test/*N=2*scan*") + glob.glob("Data/Test/NOT_CONVERGED/*N=2*scan*"),
    key=lambda s: float(s.split("beta=")[1].split("_")[0]),
)

betas, rvals, clean = [], [], []
for fn in files:
    with h5py.File(fn, "r") as f:
        beta = float(fn.split("beta=")[1].split("_")[0])
        g00 = f["results/Re_correlation/0-0"][0]
        g01 = f["results/Re_correlation/0-1"][0]
        mid = len(g00) // 2
        r = g01[mid] / g00[mid]
        negev = f["runtimedata"].attrs["negative_eigenvalue_ratio_list"][-1]
        betas.append(beta)
        rvals.append(abs(r))
        # runs whose output left the physical region (|r|>1 or big eigenvalue truncation)
        clean.append(abs(r) <= 1.05 and abs(negev) < 0.1)

betas = np.array(betas)
rvals = np.array(rvals)
clean = np.array(clean)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.6))

# --- panel 1: order parameter vs beta ---
ax1.axvspan(15, 25, color="0.9", zorder=0)
shown = np.where(clean, rvals, 1.0)
ax1.plot(betas[clean], shown[clean], "o-", color="C0", label="converged sampling")
ax1.plot(betas[~clean], shown[~clean], "s", mfc="none", color="C3", ms=9,
         zorder=5, clip_on=False, label="runaway (unphysical)")
ax1.axhline(1.0, color="0.5", lw=0.8, ls=":")
ax1.set_xlim(betas.min(), betas.max())
ax1.set_ylim(0, 1.1)
ax1.set_xlabel(r"$\beta$")
ax1.set_ylabel(r"$|g^{zz}_{01}(\beta/2)\,/\,g^{zz}_{00}(\beta/2)|$")
ax1.legend(loc="upper left")

# --- panel 2: freezing condition across cluster sizes ---
lam = np.array([3.247, 2.375, 1.125])     # pattern gains N=9, N=4, N=2
beta_c = np.array([2.2, 4.2, 18.0])       # observed transitions
beta_c_err = np.array([[0.2, 0.2, 2.0], [0.5, 0.3, 2.0]])
ax2.errorbar(beta_c, 1.0 / lam, xerr=beta_c_err, fmt="o", color="C0")
bb = np.linspace(2.0, 20.0, 100)
a_fit = np.polyfit(np.log(beta_c), np.log(1.0 / lam), 1)
ax2.plot(bb, np.exp(a_fit[1]) * bb ** a_fit[0], "-", color="C1",
         label=r"$f(\beta)\propto\beta^{%.2f}$" % a_fit[0])
for x, y, lab in zip(beta_c, 1.0 / lam, ["$N=9$", "$N=4$", "$N=2$"]):
    ax2.annotate(lab, (x, y), textcoords="offset points", xytext=(6, -10))
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlim(2.0, 20.0)
ax2.set_xlabel(r"$\beta_c$")
ax2.set_ylabel(r"$1/\lambda$")
ax2.legend(loc="upper left")

fig.tight_layout()
fig.savefig("Plots/N2_uncoupled_freezing_scan.pdf")
fig.savefig("Plots/N2_uncoupled_freezing_scan.png", dpi=200)
print("saved Plots/N2_uncoupled_freezing_scan.{pdf,png}")
print("fitted exponent a = %.3f" % a_fit[0])
