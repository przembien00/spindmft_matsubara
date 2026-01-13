import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import import_data as id

def quad(x,a,p):
    return a*x**p

# undisturbed signal:
root_folder = "selfcons"
model_frame = "dipole_lab"
symmetry_type = np.array(("B","C","D"))
phys_data = "no_pulse"
extension = np.array(("hs","l31")) # hs: home station, l31: office computer
colors = np.array(("orange","blue","k"))
lines = np.array(("-","--"))
markers = np.array(("x","v"))
ne = 1#len(extension)
ns = 3#len(symmetry_type)

tpot = np.array((0,1,2,3,4,5,6,7,8,9,10,11,12))
num_steps = np.zeros(len(tpot))
time = np.zeros(len(tpot))

data = np.zeros(((ne,ns,4)))
s = np.array(((3,3,3),(3,3,3)))

for stype in range(0,ns):
    for work_station in range(0,ne):
        for i in range(0,len(tpot)):
            param_str = id.read_PUL_params(root_folder,model_frame,symmetry_type[stype],phys_data,extension[work_station] + "CT_t_sq2t" + str(tpot[i]))
            time[i] = param_str[16].split(' ')[6] # in seconds
            num_steps[i] = param_str[2].split(' ')[0] 

        #popt,pcov = curve_fit(quad,num_steps[s[work_station][stype]:],time[s[work_station][stype]:])
        #a = popt[0]
        #p = popt[1]
        #da = np.sqrt(pcov[0,0])
        #dp = np.sqrt(pcov[1,1])
        #print("=====================================================")
        #print(extension[work_station],",",symmetry_type[stype],":")
        #print("a =",a,"+-",da)
        #print("p =",p,"+-",dp)
        #data[work_station][stype] = ([a,da,p,dp])

        #x = np.linspace(3*10**1,3*10**3,1000)
        #plt.plot(x,quad(x,a,p),ls=lines[work_station],color=colors[stype],label=r'data ' + f'{extension[work_station]}' + ', ' + f'{symmetry_type[stype]}')
        plt.plot(num_steps,time,markersize=4,markeredgewidth=0.5,marker=markers[work_station],ls='None',color=colors[stype])#,label=r'data for ' + f'{extension[work_station]}' + ', ' + f'{symmetry_type[stype]}')

#np.savetxt("../steps_data_hs.txt",data[0],header="  a  ||  da  ||  p  ||  dp")
#np.savetxt("../steps_data_l31.txt",data[1],header="  a  ||  da  ||  p  ||  dp")

plt.xscale('log')
plt.yscale('log')
#plt.xlim(3*10**1,3*10**3)
plt.xlabel(r'\# time steps')
plt.ylabel(r'$\mathrm{time} [\mathrm{s}]$')
plt.legend(loc='lower right')
plt.savefig('CT_time_steps.pdf')
