import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys
sys.path.append("../python_libs") 

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

beta_array = [0.5, 1.5, 3.5]
h_array_small = [0.001, 0.01, 0.02, 0.05]
h_array = np.arange(0.15, 0.45, 0.05)
h_array = np.concatenate((h_array_small, h_array))
# h_array = [0.001, 0.005, 0.05, 0.1, 0.15, 0.2, 0.25, 0.4, 0.45, 0.55, 0.6]
# h_vals = np.concatenate((-h_array[::-1], h_array))
for beta in beta_array:
    mags = []
    for h in h_array:
        all, disc = ImportData_spinDMFT("ISO",physical_data=f"JL=-2__beta={beta:.2g}__h=z_h_abs={h:.2g}",project="Magnetization",extension="")
        mags.append(all['results'].attrs['S_z'])
    # mags = np.concatenate(( -np.array(mags)[::-1], np.array(mags)))
    plt.plot(h_array, mags, label=rf"$\beta J_Q$={beta:.2g}")
plt.xlabel("B")
plt.ylabel("M")
plt.legend()
plt.savefig("Plots/Mag_vs_B.pdf")