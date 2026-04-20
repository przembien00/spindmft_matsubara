import numpy as np
import matplotlib.pyplot as plt
import h5py as h5



def ImportData(physical_data, project_name = ""):
    # process the inserted data:
    root_folder = "Data"

    # determine the folder:
    foldername = root_folder + "/"
    if project_name != "":
        foldername += project_name + "/"

    # determine file and return data:
    filename = foldername + physical_data + ".hdf5"
    all = h5.File( filename, 'r' )

    # discretization
    params =  all['parameters']
    disc = np.linspace(0., 1, 2*params.attrs['num_TimePoints'])

    return all, disc

def ImportData_ED(physical_data, project_name = ""):
    # process the inserted data:
    root_folder = "Data/ED_new"

    # determine the folder:
    foldername = root_folder + "/"
    if project_name != "":
        foldername += project_name + "/"

    # determine file and return data:
    filename = foldername + physical_data + ".hdf5"
    all = h5.File( filename, 'r' )

    # discretization
    params =  all['parameters']
    disc = np.linspace(0., 1, params.attrs['num_TimePoints'])

    return all, disc

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

N=12
beta_array = [0.2, 0.5, 1, 1.5, 2., 2.5]
h_z = 0.5
foldername = f"Plots/Random_h_z={h_z:.1g}_N=12/"

markers = ['v', '^', 's', 'x', 'D', 'p']
fig, (ax_xx, ax_xy, ax_zz) = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(8, 9), gridspec_kw={'hspace': 0})

i=0
for beta in beta_array:

    all, times = ImportData(f"Random_Couplings/ISO__Random__N={N}__beta={beta:.2g}__h_z={h_z:.1g}__numConfigs=250")
    # all, times = ImportData(f"ISO__Square_NN_PBC_N=20__beta={beta:.2g}__h_z={h_z:.1g}__rescale=-0.5")
    all_sdmft, times_sdmft = ImportData_spinDMFT("ISO", physical_data=f"beta={beta:.2g}__h=z_h_abs={h_z:.1g}", project="spinDMFT", extension="")
    Re_G = np.array( [ gab for gab in all['results']['Re_correlation']] )
    Im_G = np.array( [ gab for gab in all['results']['Im_correlation']] )
    Re_G_err = np.array( [ gab for gab in all['results']['Re_stddev']] )
    Im_G_err = np.array( [ gab for gab in all['results']['Im_stddev']] )
    Re_G_sdmft = np.array( [ gab for gab in all_sdmft['results']['Re_correlation']] )
    Im_G_sdmft = np.array( [ gab for gab in all_sdmft['results']['Im_correlation']] )
    G_xx = np.concatenate((Re_G[0],np.flip(Re_G[0])))
    G_xx_err = np.concatenate((Re_G_err[0],np.flip(Re_G_err[0])))
    G_zz = np.concatenate((Re_G[3],np.flip(Re_G[3])))
    G_zz_err = np.concatenate((Re_G_err[3],np.flip(Re_G_err[3])))
    G_xy = np.concatenate((Im_G[1],-np.flip(Im_G[1])))
    G_xy_err = np.concatenate((Im_G_err[1],np.flip(Im_G_err[1])))
    
    ax_xx.plot(times_sdmft, Re_G_sdmft[0], label = rf'spinDMFT, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, color='limegreen')
    ax_xx.errorbar(times, G_xx, yerr=G_xx_err, ls=':', label = rf'Chebyshev, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, color='darkviolet', errorevery=5)
    ax_xy.plot(times_sdmft, Im_G_sdmft[1], label = rf'spinDMFT, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, color='limegreen')
    ax_xy.errorbar(times, G_xy, yerr=G_xy_err, ls=':', label = rf'Chebyshev, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, color='darkviolet', errorevery=5)
    ax_zz.plot(times_sdmft, Re_G_sdmft[3], label = rf'spinDMFT, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, color='limegreen')
    ax_zz.errorbar(times, G_zz, yerr=G_zz_err, ls=':', label = rf'Chebyshev, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=25, color='darkviolet', errorevery=5)
    i+=1

marker_p = []
for m in markers:
    sup_fig, sup_ax = plt.subplots()
    picture, = sup_ax.plot([1], marker=m, color="black")
    marker_p.append(picture)
label_list = []
for beta in beta_array:
    label_list.append(rf'$\beta J_Q$={beta:.2g}')

for ax in [ax_xx, ax_xy, ax_zz]:
    ax.set_xlabel(r'$\tau$/$\beta$')
    ax.set_xlim(0, 1)
ax_zz.legend(marker_p, label_list, loc=3)


ax_xx.set_ylabel(r'$g_{xx}$($\tau$)')
ax_xy.set_ylabel(r'Im $g_{xy}$($\tau$)')
ax_zz.set_ylabel(r'$g_{zz}$($\tau$)')
# AFM:
# ax_xx.set_ylim(0.21, 0.252)
# ax_zz.set_ylim(0.21, 0.252)

# Random:
ax_xx.set_ylim(0.165, 0.252)
ax_zz.set_ylim(0.165, 0.252)
ax_xx.text(0.95, 0.17, '(a)')
ax_xy.text(0.95, -0.118, '(b)')
ax_zz.text(0.95, 0.17, '(c)')

# FM:
# ax_xx.text(0.95, 0.085, '(a)')
# ax_zz.text(0.95, 0.2318, '(c)')
# ax_xy.text(0.95, -0.225, '(b)')


fig.savefig(foldername+"stack.pdf", dpi=1000)