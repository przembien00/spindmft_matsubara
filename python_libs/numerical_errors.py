import numpy as np
import h5py as h5
import import_data as id
from plot_algorithms import DiscretizationError_pow2
from colorama import Fore, Back, Style, init

# rounding
r2 = lambda value : np.around(value,2) # round values for terminal output
r3 = lambda value : np.around(value,3) # round values for terminal output
r4 = lambda value : np.around(value,4) # round values for terminal output
r5 = lambda value : np.around(value,5) # round values for terminal output

# colors:
red, blue, green = Fore.RED, Fore.BLUE, Fore.GREEN

# return colored string
def colored(str,color):
    return color + str + Style.RESET_ALL

# return colored string if condition is fulfilled
def colorif(str,color,condition):
    if condition: 
        return colored(str,color)
    return str

# sigma : statistical error
# dQ    : time discretization error
# dI    : iteration error (only in CspinDMFT)
# epsilon : error tolerance (epsilon=1 means everything is tolerated)


# =============== spinDMFT =============== 
def max_std_spinDMFT( all, epsilon=1, printit=True ):
    M = all['parameters'].attrs['num_Samples']
    sigma_max = 4/np.sqrt(M)*np.max(all['runtimedata']['correlation_sample_std'])
    if printit:
        print("max[ σ(G^alphabeta(t)) ] =", colorif(str(sigma_max),red,sigma_max>epsilon))
    return sigma_max

def max_dQ_spinDMFT( all, all_hd, epsilon=1, printit=True ):
    G, G_hd = id.get_G_spinDMFT(all), id.get_G_spinDMFT(all_hd)
    dQ_max = 1/4 * np.max( [np.max(np.abs( DiscretizationError_pow2( Gab, Gab_hd ))) for Gab, Gab_hd in zip(G,G_hd)] )
    if printit:
        print("max_{ab,t}[ ΔQ^ab(t) ] =", colorif(str(dQ_max),red,dQ_max>epsilon))
    return dQ_max

def max_errors_spinDMFT( all, all_hd, epsilon=1, printit=True ):
    sigma_max, dQ_max = max_std_spinDMFT( all, printit=False ), max_dQ_spinDMFT( all, all_hd, printit=False )
    if printit:
        print("[σ,ΔQ] = [", colorif(str(sigma_max),red,sigma_max>epsilon), colorif(str(dQ_max),red,dQ_max>epsilon), "]")
    return sigma_max, dQ_max


# =============== CspinDMFT =============== 
def max_std_CspinDMFT( all, epsilon=1, epsilon2=1, printit=True, cindex=0 ):
    if all['parameters'].attrs['adaptive_num_SamplesPerCore'] == 1:
        M = all['parameters'].attrs['num_Cores']*all['runtimedata'].attrs['dynamical sample size per core (including next iteration)'][-2] # second to last = actual sample size
    else:
        M = all['parameters'].attrs['num_Samples']
    sample_stds = all['runtimedata']['correlation_stds'] 
    sigma_maxs = [4/np.sqrt(M) * np.max(sample_stds[key]) for key in sample_stds.keys()]
    if printit:
        cindex_entry = str(cindex)+": "+colorif(str(sigma_maxs[cindex]),red,sigma_maxs[cindex]>epsilon)
        print( "max_{ab,t}[ σ(G^ab_ij(t)) ] = [", " ,".join([cindex_entry] + [str(index)+": "+colored(str(sigma),red) for index,sigma in enumerate(sigma_maxs) if index!=cindex and sigma>epsilon2]), "]")
    return sigma_maxs

def max_dQ_CspinDMFT( all, all_hd, epsilon=1, epsilon2=1, printit=True, cindex=0 ):
    G, G_hd = id.get_G_CspinDMFT(all), id.get_G_CspinDMFT(all_hd)
    dQ_maxs = [1/4 * np.max(dQ_ij_max) for dQ_ij_max in [[np.max(np.abs( DiscretizationError_pow2( G_ab_ij, G_hd_ab_ij ))) for G_ab_ij, G_hd_ab_ij in zip(G_ij,G_hd_ij)] for G_ij,G_hd_ij in zip(G,G_hd)] ]
    if printit:
        colorstr = lambda index,dQ : str(index)+": "+colorif(str(dQ),red,dQ>epsilon2 or dQ>epsilon and index==cindex)
        print( "max_{ab,t}[ ΔQ^ab_ij(t) ] == [", " ,".join([colorstr(index,dQ) for index,dQ in enumerate(dQ_maxs)]), "]")
    return dQ_maxs

def max_dI_CspinDMFT( all, R="abs", epsilon=1, printit=True ):
    dI = 4*all['runtimedata'].attrs[R+'_iteration_error (R)'][-1]
    if printit:
        print( "max_{ab,ij}[ ΔI_%s(abij) ] ="%R, colorif(str(dI),red,dI>epsilon) )
    return dI

def max_errors_CspinDMFT( all, all_hd, R="abs", epsilon=1, epsilon2=1, printit=True, cindex=0, name="" ):
    sigma_maxs, dQ_maxs, dI = max_std_CspinDMFT( all, printit=False ), max_dQ_CspinDMFT( all, all_hd, printit=False ), max_dI_CspinDMFT( all, R=R, printit=False )
    if printit:
        sigma, dQ = sigma_maxs[cindex], dQ_maxs[cindex] # errors of important correlation
        print(name+"[ΔI,σ,ΔQ]["+str(cindex)+"] -- [σ][rest] -- [ΔQ][rest] = ["+colorif(str(r5(dI)),red,dI>epsilon),colorif(str(r5(sigma)),red,sigma>epsilon),colorif(str(r5(dQ)),red,dQ>epsilon)+"] -- ["+", ".join([str(i)+": "+colored(str(r5(dsi)),red) for i,dsi in enumerate(sigma_maxs) if dsi>epsilon2])+"] -- ["+", ".join([str(i)+": "+colored(str(r5(dQi)),red) for i,dQi in enumerate(dQ_maxs) if dQi>epsilon2])+"]") 
    return sigma_maxs, dQ_maxs, dI


# =============== nl-spinDMFT =============== 
# same as for CspinDMFT but without iteration error
def max_errors_nlspinDMFT( all, all_hd, epsilon=1, printit=True, name="" ):
    sigma_maxs, dQ_maxs = max_std_CspinDMFT( all, printit=False ), max_dQ_CspinDMFT( all, all_hd, printit=False )
    sigma, dQ = sigma_maxs[0], dQ_maxs[0] # errors of important correlation
    if printit:
        print(name+"[σ,ΔQ] =","["+colorif(str(r5(sigma)),red,sigma>epsilon),colorif(str(r5(dQ)),red,dQ>epsilon)+"]") 
    return sigma, dQ