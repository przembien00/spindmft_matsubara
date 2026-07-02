import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

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

coarse_beta_array = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
fine_beta_array = sorted(coarse_beta_array + [1.6, 1.7, 1.8, 1.9, 2.1, 2.2, 2.3, 2.4])
h_abs_extensions = {0.1: "", 0.05: "", 0.01: "chain", 0.005: "chain", 0.001: "chain", 0.0: "chain"}
fine_h_abs = {0.1, 0.01, 0.001, 0.0}
JL = -2

for h_abs, extension in h_abs_extensions.items():
    beta_array = fine_beta_array if h_abs in fine_h_abs else coarse_beta_array
    S_z = []
    S_z_err = []
    for beta in beta_array:
        all = ImportData_spinDMFT("ISO",physical_data=f"JL={JL}__beta={beta:.2g}__h=z_h_abs={h_abs:.2g}",project="Mag_MH",extension=extension)
        S_z.append(-all['results'].attrs['S_z'])
        S_z_err.append(all['runtimedata'].attrs['S_z_sample_std'])
    plt.errorbar(beta_array, S_z, yerr=S_z_err, marker='o', label=rf"$h_z$={h_abs:.2g}")

plt.xlabel(r"$\beta J_Q$")
plt.ylabel(r"$-\left<\mathbf{S}^z_0\right>$")
plt.xlim(min(coarse_beta_array), max(coarse_beta_array))
plt.legend()
plt.savefig("Plots/Mag_MH_vs_beta.pdf")
