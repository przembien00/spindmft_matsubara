import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline


def ImportData_spinDMFT( spin_model, physical_data = "", project = "", selfcons = True, extension = "" ):

    # process the inserted data:
    root_folder = "Data/"
    if physical_data != "":
        physical_data = "__" + physical_data
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = root_folder
    if not selfcons:
        foldername += "noselfcons/"
    if project != "":
        foldername += project + "/"

    # determine file and return data:
    filename = foldername + "spinmodel=" + spin_model + physical_data + extension + ".hdf5"
    all = h5.File( filename, 'r' )

    return all


JL = -2
h_abs = 0.001
beta_array = np.array([round(1.5 + 0.1*i, 2) for i in range(21)])
beta_array_dense = np.linspace(beta_array.min(), beta_array.max(), 500)

M_array = []
M_err = []
for beta in beta_array:
    all = ImportData_spinDMFT("ISO",physical_data=f"JL={JL}__beta={beta:g}__h=z_h_abs={h_abs:g}",project="Magnetization_new")
    M_array.append(all['results'].attrs['S_z'])
    M_err.append(all['runtimedata'].attrs['S_z_sample_std'])
M_array = -np.array(M_array)
M_err = np.array(M_err)

M_spl = UnivariateSpline(beta_array, M_array, s=0, k=5)
d2M_spl = M_spl.derivative(n=2)

roots = d2M_spl.roots()
beta_c = roots[np.where(roots > 2.2)][0]

print("beta_c = ", beta_c)

plt.errorbar(beta_array, M_array, yerr=M_err, fmt='o', label="Data")
plt.plot(beta_array_dense, M_spl(beta_array_dense), label="Spline", linestyle="--")
plt.axvline(beta_c, color="black", linestyle="--", label=rf"$\beta_c J_Q$={beta_c:.3g}")
plt.legend()
plt.xlabel(r"$\beta J_Q$")
plt.ylabel(r"$-\left<\mathbf{S}^z_0\right>$")
plt.xlim(beta_array.min(), beta_array.max())
plt.savefig("Plots/CT_new_vs_beta.pdf")
