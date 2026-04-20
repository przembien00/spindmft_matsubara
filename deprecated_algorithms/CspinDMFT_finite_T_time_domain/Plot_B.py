import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys
sys.path.append("../python_libs") 

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


beta_array = [19, 20, 21, 22]


markers = ['v', '^', 's', 'x', 'D', '+']


i=0
for beta in beta_array:
    all, disc = ImportData_spinDMFT("ISO",physical_data=f"JQ={beta:.2g}__beta=1",project="Iterative_Init",extension="")
    # all_c, disc_c = ImportData_spinDMFT("ISO",physical_data=f"JL=2__beta={beta:.2g}__h=z_h_abs=2",project="Complex_mfs",extension="")

    # all_cet, times = ImportData(f"ISO__Square_NN_PBC_N=20__beta={beta:.2g}__h_z=2__rescale=-0.5", project_name="Chebyshev_B_field")
    G = np.array([gab for gab in all['results']['Re_correlation']][0])
    # G_c = np.array([gab for gab in all_c['results']['Im_correlation']][1])
    # G_cet = np.array( [ gab for gab in all_cet['results']['Im_correlation']][1] )
    # g_cet = np.concatenate((-G_cet, np.flip(G_cet)))
    # C = G[0][0:int(len(G[0])/2 + 1)]
    # C = np.append(C, np.flip(C[1:]))
    plt.plot(disc,G, ls="-", label=rf'spinDMFT, $\beta J_Q$={beta:.2g}')
    # plt.plot(disc_c, G_c, ls="--", label=rf'spinDMFT Complex MF, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=20)
    # plt.plot(disc.t/beta, g_cet, ls=":", label=rf'CET N=20, $\beta J_Q$={beta:.2g}', marker=markers[i], markevery=20)
    i+=1

plt.legend()
plt.legend(fontsize=8)
plt.xlabel(r'$\tau/\beta$')
plt.ylabel(r'$g_{xx}(\tau)$')
# plt.xlim(left=0,right=1.)
# plt.ylim(bottom=0.6,top=1)
plt.savefig('Plots/Plot_22.pdf', dpi=1000)