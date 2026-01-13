import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys
sys.path.append("../python_libs") 
import scipy.fft as fft
import scipy.linalg as lin

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

def index(n):
    if n==0:
        return 0
    if n>0:
        return 2*n -1
    else:
        return -2*n

S = 0.5*np.array((np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),  np.array([[1, 0], [0, -1]])))

def compute_propagator(V_old, V_new, dt):
    T = np.array([V_new, V_old])
    result = np.eye(2, dtype=complex)
    for b in CFET:
        F = 0.5 * dt * b @ T
        U = np.cosh(lin.norm(F)*dt/2) * np.eye(2, dtype=complex)
        for i in range(3):
            U += - 2 * np.sinh(lin.norm(F)) * F[i] * S[i] / lin.norm(F)
        result = result @ U
    return result

CFET = 0.5 * np.array([[11/40+60/87+35/50, 11/40-60/87+35/50],[9/20-35/25, 9/20-35/25],[11/40-60/87+35/50, 11/40+60/87+35/50]])

beta = 0.6
h_z = 2

# all, disc = ImportData_spinDMFT("ISO",physical_data=f"JQ={beta}__beta=1",project="Iterative_Init",extension="")
all, disc = ImportData_spinDMFT("ISO",physical_data=f"beta={beta}__h=z_h_abs={h_z}",project="Mag_Field",extension="")
disc = disc[:-1] * beta
G_xx = np.array([gab for gab in all['results']['Re_correlation']][0])[:-1]
G_zz = np.array([gab for gab in all['results']['Re_correlation']][3])[:-1]
S_z = all['results'].attrs['S_z']

# EV = []
# for m, G in [(0., G_xx), (S_z, G_zz)]:
#     Eigvals = np.empty(len(disc))
#     for n in range(len(disc)):
#         integral = 0.
#         for nt1,t1 in enumerate(disc):
#                 for nt2,t2 in enumerate(disc):
#                     integral += (G[np.abs(nt1-nt2)]-m**2) * np.exp( 1j * 2 * np.pi * n * (t1 - t2) / beta )
#         Eigvals[n] = np.real(integral) * (disc[1]-disc[0])/beta
#     indices = np.argsort(Eigvals)
#     EV.append(np.sort(Eigvals))

EV = []
for m, G in [(0., G_xx), (S_z, G_zz)]:
    G_n = np.real(fft.fft(G-m**2))
    EV.append(G_n)
print(len(EV[0]))
t_arr = np.arange(len(disc))
tpoints = t_arr/len(t_arr)

O = np.empty((len(disc),len(disc)))
O[:,0] = np.ones(len(t_arr))/len(t_arr)**0.5
for n in t_arr[1:len(t_arr)//2+1]:
    O[:,n] = np.cos(2*np.pi*tpoints*n) / (len(t_arr)/2)**0.5
for n in t_arr[len(t_arr)//2+1:len(t_arr)]:
    O[:,n]  = -np.sin(2*np.pi*tpoints*n) / (len(t_arr)/2)**0.5
O_M = O

D = np.diag(EV[1])
Cov_M = O_M @ D @ O_M.T

Cov_t = np.empty((len(disc),len(disc)))
for t1 in t_arr:
    for t2 in t_arr:
        Cov_t[t1,t2] = G_xx[np.abs(t1-t2)]
EV_diag = []
eig, eigv = lin.eigh(Cov_t)
EV_diag.append(eig)
Cov_t = np.empty((len(disc),len(disc)))
for t1 in t_arr:
    for t2 in t_arr:
        Cov_t[t1,t2] = G_zz[np.abs(t1-t2)] - S_z**2
eig, eigv = lin.eigh(Cov_t)
EV_diag.append(eig)
O = eigv

Cov_t = np.empty((len(disc),len(disc)))
for t1 in t_arr:
    for t2 in t_arr:
        Cov_t[t1,t2] = G_zz[np.abs(t1-t2)] - S_z**2

print(np.max(np.abs( Cov_t - Cov_M )))

N_samples = 1000

G_xx_computed = np.zeros(len(disc), dtype=complex)

Z = 0.

# for n in range(N_samples):
#     print(n)
#     V_x = np.random.multivariate_normal(mean=np.zeros(len(EV[0])), cov=np.diag(EV[0]) )
#     V_y = np.random.multivariate_normal(mean=np.zeros(len(EV[0])), cov=np.diag(EV[0]) )
#     V_z = np.random.multivariate_normal(mean=np.zeros(len(EV[0])), cov=np.diag(EV[1]) )
#     for V in [V_x, V_y, V_z]:
#         V = [np.cos(2*np.pi*k)V[k] for k in tpoints ] 
#         # V_x = O @ V_x
#         # V_y = O @ V_y
#         # V_z = O @ V_z
#     # V_x = fft.irfft(V_x, n=len(disc))
#     # V_y = fft.irfft(V_y, n=len(disc))
#     # V_z = fft.irfft(V_z, n=len(disc))

#     shortstep_U_list = []
#     for i in range(1, len(disc)):
#         V_new= np.array([V_x[i], V_y[i], V_z[i] + h_z])
#         V_old= np.array([V_x[i-1], V_y[i-1], V_z[i-1] + h_z])
#         U = compute_propagator(V_old, V_new, disc[1]-disc[0])
#         shortstep_U_list.append(U)
#     V_new= np.array([V_x[0], V_y[0], V_z[0] + h_z])
#     V_old= np.array([V_x[-1], V_y[-1], V_z[-1] + h_z])
#     U = compute_propagator( V_old, V_new, disc[1]-disc[0] )
#     shortstep_U_list.append(U)
#     U_list = [np.eye(2)]
#     U_inv_list = [np.eye(2)]
#     for i in range(len(shortstep_U_list)):
#         U_list.append( shortstep_U_list[i] @ U_list[-1] )
#         U_inv_list.append( U_inv_list[-1] @ shortstep_U_list[-i] )
#     U_list = np.array(U_list)
#     U_inv_list = np.array(U_inv_list[::-1])
#     # print(U_list @ U_inv_list)
#     Z += np.linalg.trace(U_list[-1])
#     G_xx_sample = []
#     for i in range(len(disc)+1):
#         G_xx_sample.append( np.linalg.trace( U_inv_list[i] @ S[0] @ U_list[i] @ S[0] ) )
#     G_xx_computed += np.array(G_xx_sample, dtype=complex)
# disc = np.append(disc, beta)
# plt.plot(disc, np.real(G_xx_computed)/Z, ls="-", label=rf'Sampled, $\beta J_Q$={beta:.2g}')
# plt.plot(disc[:-1], G_xx)


# plt.plot(t_arr/(len(t_arr)), O[:,i], ls="--")
# G_iw = fft.fft(G)
# Cov_appr = O.T @ Cov_t @ O
# print(np.sort(np.diag(Cov_appr)))
# print(np.sort(Eigvals))
# print(np.sort(eig))
# for i,A in enumerate([Eigvals, eig, np.diag(Cov_appr)]):
    # plt.plot(t_arr/(len(t_arr)), np.sort(A), '.', label=f"Method {i+1}")
# plt.legend()
# print(np.max(np.abs( Cov_appr - np.diag(np.diag(Cov_appr)) )) )

# plt.plot(np.arange(len(EV_diag[0])), EV_diag[0], 'o')
# plt.plot(np.arange(len(EV[0])), EV[0], 'o')
# print(EV[0])
# i = -1
# plt.plot(disc, O_M[:,i])
# plt.plot(disc, O[:,i])




plt.savefig(f"Plots/EV.pdf")