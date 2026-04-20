import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline


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


def f(x, M, a, b, M0):
    return M * np.tanh(a * (x - b)) + M0

def w(x, Tc, k):
    return 1/(1 + np.exp(k * (1 - x/Tc)))

def g(x, Tc, alpha, A):
    M = np.where(x < Tc, np.abs(1-x/Tc)**alpha * A, 0)
    return M

def logT(x, Tc):
    return np.where(x < Tc, np.log(1-x/Tc), np.log(x))

beta_array = np.array([ 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5])
beta_array_dense = np.linspace(1.6, 3.5, 500)
M_array = []
M_array_2 = []
for beta in beta_array:
    all, disc = ImportData_spinDMFT("ISO",physical_data=f"JL=-2__beta={beta:.2g}__h=z_h_abs=0",project="Magnetization",extension="")
    M_array.append(all['results'].attrs['S_z'])


M_array = -np.array(M_array)
p0 = [0.45, 0.4, 0.5]

# popt, pcov = curve_fit(g, 1/beta_array, -M_array, p0=p0)
# print(popt)
# print("beta_c = ", 1/popt[0])
# print(g(0.3, *popt))

# dx2M = np.diff(M_array, n=1)
M_spl = UnivariateSpline(beta_array, M_array, s=0,k=5)
d2M_spl = M_spl.derivative(n=2)

roots = d2M_spl.roots()
beta_c = roots[np.where(roots>2.2)][0]

print("Tc = ", beta_c)

plt.plot(beta_array, M_array, 'o', label="Data")
plt.plot(beta_array_dense, M_spl(beta_array_dense), label="Spline", linestyle="--")
# plt.plot(beta_array_dense, d2M_spl(beta_array_dense), label="Second derivative")
# plt.plot(beta_array[:-1], dx2M, 'x', label=r"$\Delta^2 M$")
# plt.plot(temp_array_dense, g(temp_array_dense, *popt), label="Fit", linestyle="--")
plt.legend()
# plt.plot(beta_array, f(beta_array, *popt), label="Fit", linestyle="--")
plt.xlabel(r"$\beta J_Q$")
plt.ylabel(r"$\left<\mathbf{S}^z_0\right>$")
plt.xlim(1.58, 3.52)
plt.ylim(0, 0.45)
plt.savefig("Plots/M_vs_beta.pdf")
