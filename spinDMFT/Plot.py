import numpy as np 
import matplotlib.pyplot as plt
import sys
sys.path.append("../python_libs") 
import import_data as id
plt.style.use('ggplot')


beta_array = [1.0]
# all_ed = id.ImportData_ED_imagtime(f"ISO_Square_NN_PBC_N=16__ISO_Disordered_Blockwise__rescale=0.5",  root_folder='spinDMFT/Data/ExactDiagonalization')
# time = np.linspace(0, 1, all_ed['parameters'].attrs['num_TimePoints'])

markers = ['v', '^', 's', 'x', 'D', '+']


i=0
for beta in beta_array:
    all, disc = id.ImportData_spinDMFT("ISO",physical_data=f"beta={beta:.2g}XX",project="",extension="")
    # G_ed = np.array( [ 4*gab for gab in all_ed['results'][f"{beta:.2f}"]['fluctuation'] ] )
    G = np.array([4*gab for gab in all['results']['correlation']])
    # C = G[0][0:int(len(G[0])/2 + 1)]
    # C = np.append(C, np.flip(C[1:]))
    plt.plot(disc.t/beta, G[0], ls="-", label=rf'spinDMFT, $\beta J_Q$={beta:.2g}', markevery=30, marker=markers[i])
    # plt.plot(time, G_ed[0], ls="--", label=rf'ED N=16, $\beta J_Q$={beta:.2g}')
    i+=1


plt.legend()
plt.xlabel(r'$\tau/\beta$')
plt.ylabel(r'$4g^{xx}(\tau)$')
# plt.xlim(left=0,right=1.)
# plt.ylim(bottom=0.6,top=1)
plt.savefig('spindmft_test.png', dpi=1000)