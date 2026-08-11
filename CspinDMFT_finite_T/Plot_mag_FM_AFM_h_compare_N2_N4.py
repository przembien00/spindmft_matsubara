import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

def ImportData( project, config, beta, chemshift ):
    filename = f"Data/{project}/spinmodel=ISO__config={config}__beta={beta:g}__chemshift={chemshift}.hdf5"
    if beta == 1.5:
        x_filename = f"Data/{project}/spinmodel=ISO__config={config}__beta={beta:g}__chemshift={chemshift}X.hdf5"
        try:
            return h5.File( x_filename, 'r' )
        except OSError:
            pass
    return h5.File( filename, 'r' )

def NegMagnetization( f, average_over_sites, is_afm ):
    Sz = f['results/spin_expectation_values'][:, 2]
    Sz_err = f['runtimedata/spin_expectation_values_stds'][:, 2]
    if not average_over_sites:
        return -Sz[0], Sz_err[0]
    if is_afm:
        mask = np.sign(Sz) == np.sign(Sz[0])
    else:
        mask = np.ones_like(Sz, dtype=bool)
    return -Sz[mask].mean(), np.sqrt(np.sum(Sz_err[mask]**2)) / mask.sum()

models = {
    "FM": dict(project="FM_h", chemshift="uni_0.001", is_afm=False, configs={
        "N=4": "Cubic_3D_N=4_NN_J=-0.5",
    }),
    "AFM": dict(project="AFM_h", chemshift="stag_0.001", is_afm=True, configs={
        "N=4": "Cubic_3D_N=4_NN_J=0.5",
    }),
}

J_Q = np.sqrt(3)
beta_array = [1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.5, 3, 3.5, 4, 4.5, 5]

for model_name, model in models.items():
    for label, config in model["configs"].items():
        betas = []
        S_z = []
        S_z_err = []
        missing = []
        for beta in beta_array:
            try:
                f = ImportData(model["project"], config, beta, model["chemshift"])
            except OSError:
                missing.append(beta)
                continue
            betas.append(beta)
            val, err = NegMagnetization(f, average_over_sites=(label == "N=4" or model_name == "FM"), is_afm=model["is_afm"])
            S_z.append(val)
            S_z_err.append(err)

        if missing:
            print(f"{model_name} {label}: missing beta = {missing}")

        betaJ = J_Q * np.array(betas)
        plt.errorbar(betaJ, S_z, yerr=S_z_err, marker='o', label=model_name)

plt.xlabel(r"$\beta J_Q$")
plt.ylabel(r"$M$")
plt.xlim(J_Q * min(beta_array), J_Q * max(beta_array))
plt.legend()
plt.savefig("Plots/Mag_FM_AFM_h_compare_N4_vs_T.pdf")
