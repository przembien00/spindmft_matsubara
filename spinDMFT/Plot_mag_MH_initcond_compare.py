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

beta_array = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
h_abs = 0.01
JL = -2
runs = [
    ("2x", "burn-in=500, cold start"),
    ("noburnin", "burn-in=0, cold start"),
    ("noburnin_initH01", "burn-in=0, warm start from $h_z$=0.1"),
    ("initH01", "burn-in=500, warm start from $h_z$=0.1"),
]

for extension, label in runs:
    S_z = []
    S_z_err = []
    for beta in beta_array:
        all = ImportData_spinDMFT("ISO",physical_data=f"JL={JL}__beta={beta:.2g}__h=z_h_abs={h_abs:.2g}",project="Mag_MH",extension=extension)
        S_z.append(-all['results'].attrs['S_z'])
        S_z_err.append(all['runtimedata'].attrs['S_z_sample_std'])
    plt.errorbar(beta_array, S_z, yerr=S_z_err, marker='o', label=label)

plt.xlabel(r"$\beta J_Q$")
plt.ylabel(r"$-\left<\mathbf{S}^z_0\right>$")
plt.xlim(min(beta_array), max(beta_array))
plt.legend()
plt.savefig("Plots/Mag_MH_initcond_compare.pdf")
