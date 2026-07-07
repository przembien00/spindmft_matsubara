import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

def ImportData( config, beta, project = "FM_h", chemshift = "uni_0.1" ):
    filename = f"Data/{project}/spinmodel=ISO__config={config}__beta={beta:g}__chemshift={chemshift}.hdf5"
    return h5.File( filename, 'r' )

configs = {
    "N=2": "Cubic_3D_N=2_NN_J=-0.5",
    "N=4": "Cubic_3D_N=4_NN_J=-0.5",
}
beta_array = [0.2, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]

for label, config in configs.items():
    betas = []
    S_z = []
    S_z_err = []
    missing = []
    for beta in beta_array:
        try:
            f = ImportData(config, beta)
        except OSError:
            missing.append(beta)
            continue
        betas.append(beta)
        S_z.append(f['results/spin_expectation_values'][0, 2])
        S_z_err.append(f['runtimedata/spin_expectation_values_stds'][0, 2])

    if missing:
        print(f"{label}: missing beta = {missing}")

    plt.errorbar(betas, S_z, yerr=S_z_err, marker='o', label=label)

plt.xlabel(r"$\beta J_Q$")
plt.ylabel(r"$\left<S^z_0\right>$")
plt.xlim(min(beta_array), max(beta_array))
plt.legend()
plt.savefig("Plots/Mag_FM_h0.1_N2_N4_site0_vs_beta.pdf")
