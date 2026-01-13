import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import functions

def lin(x,m,p):
    return m*x**p

# undisturbed signal:
root_folder = "selfcons_PUL"
model_frame = "dipole_lab"
symmetry_type = "B"
phys_data = "0.00"
extension = np.array(("hs","l31")) # hs: home station, l31: office computer
colors = np.array(("orange","blue"))
ne = len(extension)

Mpot = np.array((1,2,3,4,5,6))
M = np.zeros(len(Mpot))
time = np.zeros(len(Mpot))

data = np.zeros((ne,4))
s = np.array((1,3))

for work_station in range(0,ne):
    for i in range(0,len(M)):
        param_str = functions.param_file_to_strlist(root_folder,model_frame,symmetry_type,phys_data,extension[work_station] + "CT_s_10t" + str(Mpot[i]))
        time[i] = param_str[16].split(' ')[6] # in seconds
        M[i] = 10**Mpot[i]

    popt,pcov = curve_fit(lin,M[s[work_station]:],time[s[work_station]:])

    m = popt[0]
    dm = np.sqrt(pcov[0,0])
    p = popt[1]
    dp = np.sqrt(pcov[1,1])
    print("=====================================================")
    print(extension[work_station],":")
    print("m =",m,"+-",dm)
    print("p =",p,"+-",dp)

    data[work_station] = ([m,dm,p,dp])

    x = np.linspace(4,2*10**6,1000)
    plt.plot(M,time,marker='x',ls='None',color=colors[work_station],label=r'data for %s'%extension[work_station])
    plt.plot(x,lin(x,m,p),ls='-',color=colors[work_station],label=r'linear fit for %s'%extension[work_station])

np.savetxt("../sample_data.txt",data,header="  m  ||  dm  ||  p  ||  dp")

plt.xlim(4,2*10**6)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'\# samples')
plt.ylabel(r'$\mathrm{time} [\mathrm{s}]$')
plt.legend(loc='best')
plt.savefig('CT_sample_size.pdf')
