import numpy as np
import h5py as h5
import matplotlib.pyplot as plt


def ImportData_spinDMFT( spin_model, physical_data = "", project = "", selfcons = True, extension = "", postfix=None ):

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

    # discretization
    params =  all['parameters']
    disc = np.linspace(0., 1, params.attrs['num_TimePoints'])
    
    return all, disc


beta_array = [1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.3, 2.4, 2.5]

M_array = []
for beta in beta_array:
    all, disc = ImportData_spinDMFT("ISO",physical_data=f"JL=-2__beta={beta:.2g}__h=z_h_abs=0.05",project="Magnetization",extension="")
    M_array.append(all['results'].attrs['S_z'])

plt.plot(beta_array, M_array)
plt.xlabel(r"$\beta J_Q$")
plt.ylabel(r"$\left<\mathbf{S}^z_0\right>$")
# plt.xlim(1.5, 2.5)
plt.savefig("Plots/M_vs_beta.pdf")
