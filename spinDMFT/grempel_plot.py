import csv
import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py as h5
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


betas = [1, 2, 3, 5, 7, 12]

data_X = {}
data_Y = {}


for beta in betas:
    with open('Data/grempel_data.csv', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        X = []
        Y = []
        for row in reader:
            if row[f'bJ={beta} X'] != '':
                X.append( float(row[f'bJ={beta} X'] ) )
            if row[f'bJ={beta} Y'] != '':
                Y.append( float(row[f'bJ={beta} Y'] ) )
        data_X[f'beta={beta}'] = np.array(X)
        data_Y[f'beta={beta}'] = np.array(Y)

fig, ax = plt.subplots()

for beta in betas:
    if beta == 30:
        all, disc = ImportData_spinDMFT("ISO",physical_data=f"JQ=25__beta=1",project="Iterative_Init",extension="")
    else:
        all, disc = ImportData_spinDMFT("ISO",physical_data=f"JQ={beta:.2g}__beta=1",project="Iterative_Init",extension="")
    G = np.array([gab for gab in all['results']['Re_correlation']])
    ax.plot(disc, G[0], ls="-", label=rf'spinDMFT, $\beta J$={beta:.2g}', color='limegreen')
    ax.plot(data_X[f'beta={beta}'][::], data_Y[f'beta={beta}'][::], 'x', label=rf'Grempel, $\beta J$={beta}', color='darkslategrey')
    N2 = len(G[0])//2
    ax.text(0.45, G[0][N2]+0.005, rf'$\beta J_Q$={beta:.2g}')



ax.set_xlabel(r'$\tau/\beta$')
ax.set_ylabel(r'$g^{xx}(\tau)$')
# plt.legend(prop={'size': 8})
fig.savefig('Plots/Plot_grempel.pdf', dpi=1000)
    