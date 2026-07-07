import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

def ImportData( config, beta, project = "FM", chemshift = "uni_0.001" ):
    filename = f"Data/{project}/spinmodel=ISO__config={config}__beta={beta:g}__chemshift={chemshift}.hdf5"
    return h5.File( filename, 'r' )

config = "Cubic_3D_N=4_NN_J=-0.5"
beta_array = [0.2, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]

S_z = []
S_z_err = []
for beta in beta_array:
    f = ImportData(config, beta)
    S_z.append(f['results/spin_expectation_values'][0, 2])
    S_z_err.append(f['runtimedata/spin_expectation_values_stds'][0, 2])

plt.errorbar(beta_array, S_z, yerr=S_z_err, marker='o')
plt.xlabel(r"$\beta J_Q$")
plt.ylabel(r"$\left<S^z_0\right>$")
plt.xlim(min(beta_array), max(beta_array))
plt.savefig("Plots/Mag_FM_N4_site0_vs_beta.pdf")
