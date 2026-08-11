import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

def ImportData( config, beta, project = "FM_h", chemshift = "uni_0.1" ):
    filename = f"Data/{project}/spinmodel=ISO__config={config}__beta={beta:g}__chemshift={chemshift}.hdf5"
    if beta == 1.5:
        x_filename = f"Data/{project}/spinmodel=ISO__config={config}__beta={beta:g}__chemshift={chemshift}X.hdf5"
        try:
            return h5.File( x_filename, 'r' )
        except OSError:
            pass
    return h5.File( filename, 'r' )

def NegMagnetization( f, average_over_sites ):
    Sz = f['results/spin_expectation_values'][:, 2]
    Sz_err = f['runtimedata/spin_expectation_values_stds'][:, 2]
    if not average_over_sites:
        return -Sz[0], Sz_err[0]
    return -Sz.mean(), np.sqrt(np.sum(Sz_err**2)) / len(Sz)

configs = {
    "N=2": "Cubic_3D_N=2_NN_J=-0.5",
    "N=4": "Cubic_3D_N=4_NN_J=-0.5",
}
J_Q = np.sqrt(3)
beta_array = [1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.5, 3, 3.5, 4, 4.5, 5]

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
        val, err = NegMagnetization(f, average_over_sites=True)
        S_z.append(val)
        S_z_err.append(err)

    if missing:
        print(f"{label}: missing beta = {missing}")

    betaJ = J_Q * np.array(betas)
    plt.errorbar(betaJ, S_z, yerr=S_z_err, marker='o', label=label)

plt.xlabel(r"$\beta J_Q$")
plt.ylabel(r"$M$")
plt.xlim(J_Q * min(beta_array), J_Q * max(beta_array))
plt.legend()
plt.savefig("Plots/Mag_FM_h0.1_N2_N4_site0_vs_beta.pdf")
