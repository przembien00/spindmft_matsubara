import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
import os.path
from os.path import exists
import warnings
import h5py as h5

# =============== importing the data: ===============
# root_folder   : "selfcons", "expval,spincorr", "expval,measureT", ...
# model_frame   : "isotropic_lab", "dipole_lab", ...
# symmetry_type : "A", "B", "C", "D"
# phys_data     : string with the foldername containing the physical data
# extension     : any name extensions of the data folder

# Note that the matrix returns are transposed, i.e., t1 denotes the column and t2 denotes the row
# ===================================================

def ImportData_ED_imagtime( physical_data, project = "", extension = "", root_folder = "ExactDiagonalization/Data"):

    # process the inserted data:
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation() + root_folder + "/"
    if project != "":
        foldername += project + "/"

    # determine file and return data:
    filename = foldername + physical_data + extension + ".hdf5"
    all = h5.File( filename, 'r' )

    # discretization
    params =  all['parameters']
    
    return all

# =============== 1) import TTI data ===============
# algorithm with time translation invariance => vector returns
def ImportData_TTI(root_folder, model_frame, symmetry_type, phys_data, extension):

    # process the inserted data:
    root_folder = np.asarray(root_folder.split(","))
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation()
    if root_folder[0]=="selfcons":
        foldername += "DMFT_selfcons_TTI/Data/" + model_frame + "/Type" + symmetry_type + "/" + phys_data + extension
    elif root_folder[0]=="expval":
        foldername += "DMFT_expval_TTI/Data/" + root_folder[1] + "/" + model_frame + "/Type" + symmetry_type + "/" + phys_data + extension
    elif root_folder[0]=="noise":
        foldername += "Gaussian_noise/Data/Type" + symmetry_type + "/" + phys_data + extension

    # return the data:
    return Read_TTI(foldername,symmetry_type), ReturnTime( foldername + "/params.txt" )


# =============== 2) import PUL data ===============
# algorithm without time translation invariance => matrix or vector returns 
def ImportData_PUL(root_folder, model_frame, symmetry_type, phys_data, extension):

    # process the inserted data:
    root_folder = np.asarray(root_folder.split(","))
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation()
    if root_folder[0]=="selfcons":
        foldername += "DMFT_selfcons_PUL/Data/" + model_frame + "/Type" + symmetry_type + "/" + phys_data + extension
    elif root_folder[0]=="expval":
        foldername += "DMFT_expval_PUL/Data/" + root_folder[1] + "/" + model_frame + "/Type" + symmetry_type + "/" + phys_data + extension

    # return the data:
    if root_folder[0]=="expval":
        if root_folder[1]=="measureT":
            return ReadVector_PUL(foldername), ReturnTime( foldername + "/params.txt" )[-1]    
        else:
            return ReadTwotime_PUL(foldername,symmetry_type), ReturnTime( foldername + "/params.txt" )
    else: 
        return ReadTwotime_PUL(foldername,symmetry_type), ReturnTime( foldername + "/params.txt" )


# =============== 3) import Dimer data ===============
# algorithm with time translation invariance => vector returns
def ImportData_Dimer(model_frame, symmetry_type, phys_data, extension):

    # process the inserted data:
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation() + "Dimer_DMFT_selfcons_TTI/Data/" + model_frame + "/Type" + symmetry_type + "/" + phys_data + extension

    # return the data:
    return Read_Dimer(foldername,symmetry_type), ReturnTime( foldername + "/params.txt" )


# =============== 4) import Nmer data ===============
# algorithm with time translation invariance => vector returns
def ImportData_Nmer(model_frame, symmetry_type, phys_data, N, which, extension):

    # process the inserted data:
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation() + "Dimer_DMFT_selfcons_TTI/Data/" + model_frame + "/Type" + symmetry_type + "/" + phys_data + extension

    # return the data:
    return Read_Nmer(foldername,symmetry_type,N,which), ReturnTime( foldername + "/params.txt" )


# =============== 5a) import Cluster data ===============
# algorithm with time translation invariance => vector returns
def ImportData_Cluster(model_frame, coupling_and_symmetry, phys_data, which, extension, warn_ETV = False ):

    # process the inserted data:
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation() + "Data_in_deprecated_formats/CspinDMFT/Data/" + model_frame + "/" + coupling_and_symmetry + "/" + phys_data + extension
    
    # check whether the parameters file exists and if required add ETV (eigenvalue threshold violated):
    ETV = False
    paramsfilename = foldername + "/params.txt"
    if not exists( paramsfilename ): # is the negative eigenvalues threshold is violated in the data (then the foldername obtains an extension)
        foldername = foldername + "_ETV"
        ETV = True
        paramsfilename = foldername + "/params.txt"
        if not exists( paramsfilename ):
            raise Exception( "Could not find the file : %s" %paramsfilename )

    # warn if warner is active:
    violated = False
    if ETV and warn_ETV: 
        warnings.warn( "Eigenvalue threshold violated in data set: %s" %phys_data )
        violated = True

    # return the data:
    filename = foldername + "/results_S" + str(which[0]+1) + "S" + str(which[1]+1) + ".txt"
    return Read_CorrelationTensor(filename, coupling_and_symmetry[-1]), ReturnTime( paramsfilename ), violated


# =============== 5b) import Cluster data: negative eigenvalue ratio ===============
# algorithm with time translation invariance => vector returns
def GetEigenvalueViolation_Cluster(model_frame, coupling_and_symmetry, phys_data, N, which, extension):

    # process the inserted data:
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    filename = Orientation() + "Data_in_deprecated_formats/CspinDMFT/Data/" + model_frame + "/" + coupling_and_symmetry + "/" + phys_data + extension + "/params.txt"

    file = open('%s'%filename, 'r')
    lines = file.readlines()
    for l in lines:
        line = l.strip() # ignore newline character
        if line.find( 'sum of negative EV ratios' ) != -1: # finding the right line
            value = float( line[line.find( ':' )+2:] )
            break

    # return the data:
    return value


# =============== 5c) import Cluster data: sample stds ===============
# algorithm with time translation invariance => vector returns
def ImportData_ClusterSampleSTD(model_frame, coupling_and_symmetry, phys_data, which, extension, warn_ETV = False):

    # process the inserted data:
    if extension != "":
        extension = "_" + extension

    # determine the folder and filename:
    foldername = Orientation() + "Data_in_deprecated_formats/CspinDMFT/Data/" + model_frame + "/" + coupling_and_symmetry + "/" + phys_data + extension
    
    # check whether the parameters file exists:
    ETV = False
    paramsfilename = foldername + "/params.txt"
    if not exists( paramsfilename ): # is the negative eigenvalues threshold is violated in the data (then the foldername obtains an extension)
        paramsfilename = foldername + "_ETV/params.txt"
        ETV = True
        if not exists( paramsfilename ):
            raise Exception( "Could not find the file : %s" %paramsfilename )

    # warn if warner is active:
    violated = False
    if ETV and warn_ETV: 
        warnings.warn( "Eigenvalue threshold violated in data set: %s" %phys_data )
        violated = True

    # return the data:
    filename = foldername + "/sample_stds_S" + str(which[0]+1) + "S" + str(which[1]+1) + ".txt"
    return Read_CorrelationTensor(filename, coupling_and_symmetry[-1]), ReturnTime( foldername + "/params.txt" ), violated
    

# =============== 5) import Cluster data full ===============
# algorithm with time translation invariance => vector returns
def ImportData_Cluster_full(model_frame, symmetry_type, phys_data, N, extension):

    # process the inserted data:
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation() + "Cluster_DMFT_selfcons_TTI/Data/" + model_frame + "/Type" + symmetry_type + "/" + phys_data + extension

    # return the data:
    i1 = 1
    i2 = 1
    first_correlation_set = Read_Nmer( foldername, symmetry_type, N, str(i1)+str(i2) )
    DataSet = np.zeros(( int((N*(N+1))/2), np.shape(first_correlation_set)[0], np.shape(first_correlation_set)[1] ))
    DataSet[0] = first_correlation_set
    counter = 1
    while True:
        if i2 < N :
            i2+=1
        else:
            i1+=1
            i2=i1
        if i1 > N:
            break
        double_index = str(i1)+str(i2)
        DataSet[counter] = Read_Nmer( foldername, symmetry_type, N, double_index )
        counter+=1
    
    return DataSet, ReturnTime( foldername + "/params.txt" )


# =============== 6) import Exact Diagonalization data ===============
# algorithm with time translation invariance => vector returns
def ImportData_ED_Old(model_frame, symmetry_type, phys_data, N, which, extension):

    # process the inserted data:
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation() + "Exact_Diagonalization/Data/" + model_frame + "/Type" + symmetry_type + "/" + phys_data + extension

    # return the data:
    return Read_Nmer(foldername,symmetry_type,N,which), ReturnTime_ED(foldername)


# =============== make array ===============
# this function interprets a physical data string with multiple value options generating an array of strings according to
# "text1[a,b]text2"  ===> (("text1atext2"),("text1btext2")) 
# if no separable value string is found the inserted string is returned
def MakeArray(phys_data):
    if len(phys_data[0]) == 1: # single string
        f = phys_data.find("[")
        t = phys_data.find("]")
        if f!=-1 and t!=-1: # else no separable value string in inserted array
            params = np.asarray(phys_data[f+1:t].split(','))
            output = ["" for x in range(len(params))]
            for l in range(0,len(output)):
                output[l] = phys_data[:f] + params[l] + phys_data[t+1:]
            phys_data = output
    else: # multiple strings
        output = np.array((MakeArray(phys_data[0]), MakeArray(phys_data[1])))
        for i in range(2,len(phys_data)):
            output = np.append(output,MakeArray(phys_data[i]))
        phys_data = output
    return phys_data
    
        

# =============== return plot number ===============
# returns the position of correlation "corr" in the return array A
def ReturnPlotNumber(corr):
    if corr == "xx":
        return 0,0,0,0
    elif corr == "xy":
        return np.nan,np.nan,1,1
    elif corr == "xz":
        return np.nan,np.nan,np.nan,2
    elif corr == "yy":
        return 0,0,0,3
    elif corr == "yz":
        return np.nan,np.nan,np.nan,4
    elif corr == "zz":
        return 0,1,2,5
    else:
        raise Exception("correlation %s unknown"%corr)



# =============== sub functions ================

# orientation in the repository folder tree
def Orientation( postfix = None ): 
    tree_climb = ""
    i=0 # preventing an infinite loop
    while not os.path.isfile("./%s.orientation"%tree_climb) and i<10: # climb the folder tree until the root folder is reached
        tree_climb = tree_climb + "../"
        i+=1
    if os.path.isfile("./%s.documents"%tree_climb):
        tree_climb += "../PhD_algorithms/"
    if postfix != None:
        tree_climb += postfix + "/"
    return tree_climb


# returning the linspace
def ReturnTime(filename):
    steps, Delta_t, M = np.genfromtxt( filename, unpack = 'True')
    return np.linspace(0, steps*Delta_t, int(steps)+1)

# returning the linspace
def ReturnTime_ED(foldername):
    steps, Delta_t = np.genfromtxt("%s/params.txt"%foldername, unpack = 'True')
    return np.linspace(0, steps*Delta_t, int(steps)+1)


# reading expectation values for TTI
def Read_TTI(foldername, symmetry_type):
    if symmetry_type == "A":
        gxx = np.genfromtxt("%s/results.txt"%foldername, unpack = 'True')
        return np.array((gxx))
    elif symmetry_type == "B":
        gxx,gzz = np.genfromtxt("%s/results.txt"%foldername, unpack = 'True')
        return np.array((gxx,gzz))
    elif symmetry_type == "C":
        gxx,gxy,gzz = np.genfromtxt("%s/results.txt"%foldername, unpack = 'True')
        return np.array((gxx,gxy,gzz))
    elif symmetry_type == "D":
        gxx,gxy,gxz,gyx,gyy,gyz,gzx,gzy,gzz = np.genfromtxt("%s/results.txt"%foldername, unpack = 'True')
        return np.array((gxx,gxy,gxz,gyx,gyy,gyz,gzx,gzy,gzz))


# reading expectation values for TTI
def Read_Dimer(foldername, symmetry_type):
    if symmetry_type == "B":
        gxx11,gzz11,gxx12,gzz12 = np.genfromtxt("%s/results.txt"%foldername, unpack = 'True')
        return np.array((gxx11,gzz11,gxx12,gzz12))


# reading expectation values for TTI
def Read_Nmer(foldername, symmetry_type, N, which):
    A = np.nan
    i = which[0]
    j = which[1]
    if j<i: # swap indices
        k = j 
        j = i
        i = k
    if symmetry_type == "A":
        filename = foldername + "/results_S" + str(i+1) + "S" + str(j+1) + ".txt" # in the data spins are still counted starting from 1
        gxxij = np.genfromtxt("%s"%filename, unpack = 'True')
        A = np.array(( gxxij ))
    elif symmetry_type == "B":
        filename = foldername + "/results_S" + str(i+1) + "S" + str(j+1) + ".txt" # in the data spins are still counted starting from 1
        gxxij, gzzij = np.genfromtxt("%s"%filename, unpack = 'True')
        A = np.array(( gxxij, gzzij ))
    if np.shape(A)==():
        raise Exception("data not found")
    else:
        return A

# reading a correlation tensor given a filename and the symmetry type
def Read_CorrelationTensor( filename, symmetry_type ):
    if symmetry_type == "A":
        gxx = np.genfromtxt( filename, unpack = 'True')
        A = np.array(( gxx ))
    elif symmetry_type == "B":
        gxx, gzz = np.genfromtxt( filename, unpack = 'True')
        A = np.array(( gxx, gzz ))
    elif symmetry_type == "C":
        gxx,gxy,gzz = np.genfromtxt( filename, unpack = 'True')
        A = np.array(( gxx, gxy, gzz ))
    elif symmetry_type == "D":
        gxx, gxy, gxz, gyx, gyy, gyz, gzx, gzy, gzz = np.genfromtxt( filename, unpack = 'True')
        A = np.array(( gxx, gxy, gxz, gyx, gyy, gyz, gzx, gzy, gzz ))
    else: 
        raise Exception("symmetry type %s undefined" %symmetry_type)
    return A


# reading two-time expectation values for PUL
def ReadTwotime_PUL(foldername, symmetry_type):
    if symmetry_type == "A":
        gxx = np.genfromtxt("%s/results_gxx.txt"%foldername, unpack='True')
        return np.array((gxx))
    elif symmetry_type == "B":
        gxx = np.genfromtxt("%s/results_gxx.txt"%foldername, unpack='True')
        gzz = np.genfromtxt("%s/results_gzz.txt"%foldername, unpack='True')
        return np.array((gxx,gzz))
    elif symmetry_type == "C":
        gxx = np.genfromtxt("%s/results_gxx.txt"%foldername, unpack = 'True')
        gxy = np.genfromtxt("%s/results_gxy.txt"%foldername, unpack = 'True')
        gzz = np.genfromtxt("%s/results_gzz.txt"%foldername, unpack = 'True')
        return np.array((gxx,gxy,gzz))
    elif symmetry_type == "D":
        gxx = np.genfromtxt("%s/results_gxx.txt"%foldername, unpack = 'True')
        gxy = np.genfromtxt("%s/results_gxy.txt"%foldername, unpack = 'True')
        gxz = np.genfromtxt("%s/results_gxz.txt"%foldername, unpack = 'True')
        gyy = np.genfromtxt("%s/results_gyy.txt"%foldername, unpack = 'True')
        gyz = np.genfromtxt("%s/results_gyz.txt"%foldername, unpack = 'True')
        gzz = np.genfromtxt("%s/results_gzz.txt"%foldername, unpack = 'True')
        return np.array((gxx,gxy,gxz,gyy,gyz,gzz))

# reading vector-like expectation values for PUL
def ReadVector_PUL(foldername):
    g = np.genfromtxt("%s/results.txt"%foldername, unpack='True')
    return g



# =============== 7a) Import Spin Diffusion Data ===============
# everything is saved in a single hdf5 file
def ImportData_Diff( sc_model, physical_data, project = "", extension = ""):

    # process the inserted data:
    root_folder = "Diffusion_spinDMFT/Data/"
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation() + root_folder + sc_model + "/"
    if project != "":
        foldername += project + "/"

    # determine file and return data:
    filename = foldername + physical_data + extension + ".hdf5"
    all = h5.File( filename, 'r' )

    # discretization
    params =  all['parameters']
    disc = DiscretizationCompact( params.attrs['num_PositionsPerDirection'], params.attrs['delta_r'], params.attrs['num_TimeSteps'], params.attrs['delta_t'] )

    return all, disc 

# compact storage of discretization parameters
class DiscretizationCompact:
    def __init__( self, num_pos_in_dir, delta_r, num_time_steps, delta_t ):
        self.num_pos_in_dir = num_pos_in_dir
        self.dr = delta_r
        self.num_time_steps = num_time_steps
        self.dt = delta_t
        self.r = np.linspace( 0, self.num_pos_in_dir*self.dr, self.num_pos_in_dir )
        self.t = np.linspace( 0, self.num_time_steps*self.dt, self.num_time_steps+1 )


# =============== 7b) Import Magnetization DGL Data ===============
# everything is saved in a single hdf5 file
def ImportData_MagnDGL( physical_data, project = "", extension = ""):

    # process the inserted data:
    root_folder = "Magnetization_DGL_Solver/Data"
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation() + root_folder + "/"
    if project != "":
        foldername += project + "/"

    # determine file and return data:
    filename = foldername + physical_data + extension + ".hdf5"
    all = h5.File( filename, 'r' )

    return all


# =============== 8) Import Exact Diagonalization Data ===============
# everything is saved in a single hdf5 file
def ImportData_ED( physical_data, project = "", extension = ""):

    # process the inserted data:
    root_folder = "ED_Spin_Simulation/Data"
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation() + root_folder + "/"
    if project != "":
        foldername += project + "/"

    # determine file and return data:
    filename = foldername + physical_data + extension + ".hdf5"
    all = h5.File( filename, 'r' )

    # discretization
    params =  all['parameters']
    disc = TimeDiscretizationCompact( params.attrs['num_TimeSteps'], params.attrs['delta_t'] )

    return all, disc 

# compact storage of discretization parameters
class TimeDiscretizationCompact:
    def __init__( self, num_steps, dt ):
        self.num_steps = num_steps
        self.dt = dt
        self.t = np.linspace( 0, self.num_steps*self.dt, self.num_steps+1 )

    def stretch( self, factor ):
        self.dt = factor * self.dt
        self.t = factor * self.t 

# =============== 9) Import spinDMFT Data =============== (refactored code)
# everything is saved in a single hdf5 file
def ImportData_spinDMFT( spin_model, physical_data = "", project = "", selfcons = True, extension = "", postfix=None ):

    # process the inserted data:
    root_folder = "spinDMFT/Data/"
    if physical_data != "":
        physical_data = "__" + physical_data
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation(postfix) + root_folder
    if not selfcons:
        foldername += "noselfcons/"
    if project != "":
        foldername += project + "/"

    # determine file and return data:
    filename = foldername + "spinmodel=" + spin_model + physical_data + extension + ".hdf5"
    all = h5.File( filename, 'r' )

    # discretization
    disc = TimeDiscretizationCompact( all['parameters'].attrs['num_TimeSteps'], all['parameters'].attrs['delta_t'] )

    return all, disc 

def get_G_spinDMFT( all ):
    return [4*gab for gab in all['results']['correlation']]


# =============== 9) Import spinDMFT spin 1 Data ===============
# everything is saved in a single hdf5 file
def ImportData_spinDMFT_spin1( spin_model, physical_data = "", project = "", selfcons = True, extension = ""):

    # process the inserted data:
    root_folder = "spinDMFT_spin1/Data/"
    if physical_data != "":
        physical_data = "__" + physical_data
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation() + root_folder
    if not selfcons:
        foldername += "noselfcons/"
    if project != "":
        foldername += project + "/"

    # determine file and return data:
    filename = foldername + "spinmodel=" + spin_model + physical_data + extension + ".hdf5"
    all = h5.File( filename, 'r' )

    # discretization
    params =  all['parameters']
    disc = TimeDiscretizationCompact( params.attrs['num_TimeSteps'], params.attrs['delta_t'] )

    return all, disc 

# =============== 9) Import spin echos Data ===============
# everything is saved in a single hdf5 file
def ImportData_spin_echos( spin_model, spin_echo, physical_data = "", project = "", extension = "", postfix=None):

    # process the inserted data:
    root_folder = "spin_echos/Data/"
    if physical_data != "":
        physical_data = "__" + physical_data
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation(postfix) + root_folder
    if project != "":
        foldername += project + "/"

    # determine file and return data:
    filename = foldername + "spinmodel=" + spin_model + "__echo=" + spin_echo + physical_data + extension + ".hdf5"
    all = h5.File( filename, 'r' )

    # discretization
    params =  all['parameters']
    disc = TimeDiscretizationCompact( params.attrs['num_EchoTimeSteps'], params.attrs['delta_tau'] )

    return all, disc 


# =============== 11) Import CspinDMFT Data =============== (refactored code)
# everything is saved in a single hdf5 file
def ImportData_CspinDMFT_( sc_model, spin_model, physical_data, project = "", selfcons = True, extension = "", warn_ETV = False, postfix = None ):
    # process the inserted data:
    root_folder = "CspinDMFT/Data/"
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation(postfix) + root_folder
    if not selfcons:
        foldername += "noselfcons/"
    if project != "":
        foldername += project + "/"
    if sc_model != "":
        foldername += sc_model + "/"

    # determine file and return data:
    filename = foldername + spin_model + "/" + physical_data + extension

    # check whether the parameters file exists and if required add ETV (eigenvalue threshold violated):
    ETV = False
    if not exists( filename + ".hdf5" ): # is the negative eigenvalues threshold is violated in the data (then the foldername obtains an extension)
        filename = filename + "_ETV"
        ETV = True
        if not exists( filename + ".hdf5" ):
            raise Exception( "Could not find the file : %s" %(filename + ".hdf5") )
    filename += ".hdf5"
    all = h5.File( filename, 'r' )

    # warn if warner is active:
    if ETV and warn_ETV: 
        warnings.warn( "Eigenvalue threshold violated in data set: %s" %physical_data )

    # discretization
    params =  all['parameters']
    disc = TimeDiscretizationCompact( params.attrs['num_TimeSteps'], params.attrs['delta_t'] )

    return all, disc 

def get_G_CspinDMFT( all ):
    return [[4*gab for gab in all['results']['correlation'][key]] for key in all['results']['correlation'].keys()]

# =============== 11) Import Ensemble Data from CspinDMFT_clusterization =============== 
def ImportEnsemble_CspinDMFT(subfoldername,NGamma,filename_pre=""):
    # determine the filename and import:
    filename = Orientation() + "CspinDMFT_clusterization/ClusterEnsembles/"
    if subfoldername != "":
        filename += subfoldername + "/" 
    if filename_pre != "":
        filename += filename_pre + "_"
    if NGamma == "guess":
        for NG in range(1,21):
            if os.path.isfile(filename + "NGamma=" + str(NG) + ".dat"):
                break 
        if NG > 19:
            raise RuntimeError("could not guess NGamma!")
        NGamma = NG
    filename += "NGamma=" + str(NGamma) + ".dat"
    positions = np.genfromtxt(filename, unpack='True')
    positions = [np.array([pos_in_dir[i] for pos_in_dir in positions]) for i in range(np.shape(positions)[1])] # rearrange

    # separate into cluster and environment and return 
    return positions[:NGamma], positions[NGamma:]



# =============== 12) Import nl-spinDMFT Data ===============
# everything is saved in a single hdf5 file
def ImportData_nlspinDMFT( spin_model, config, other_physical_data = None, project = None, extension = None, warn_ETV = True , postfix = None ):
    # determine the folder:
    root_folder = "nl-spinDMFT/Data/"
    foldername = Orientation(postfix) + root_folder
    if project != None:
        foldername += project + "/"

    # interpret input:
    if isinstance(spin_model,str):
        spin_model = "spinmodel=" + spin_model
    else: # spin_model is a list of two elements
        spin_model = "spinmodel=" + spin_model[0] + "__mfspinmodel=" + spin_model[1]
    if other_physical_data != None:
        if not isinstance(other_physical_data,str): # other_physical_data is a list of elements
            other_physical_data = "__".join(filter(None, other_physical_data))
        if other_physical_data != "":
            other_physical_data = "__" + other_physical_data
    else:
        other_physical_data = ""
    if extension != None:
        extension = "_" + extension
    else:
        extension = ""

    # create file str:
    filename = foldername + spin_model + "__config=" + config + other_physical_data + extension

    # check whether the parameters file exists and if required add ETV (eigenvalue threshold violated):
    ETV = False
    if not exists( filename + ".hdf5" ): # is the negative eigenvalues threshold is violated in the data (then the foldername obtains an extension)
        filename = filename + "_ETV"
        ETV = True
        if not exists( filename + ".hdf5" ):
            raise Exception( "Could not find the file : %s" %(filename + ".hdf5") )
    filename += ".hdf5"
    all = h5.File( filename, 'r' )

    # warn if warner is active:
    if ETV and warn_ETV: 
        warnings.warn( "Eigenvalue threshold violated in data set: %s" %filename )

    # discretization
    params =  all['parameters']
    disc = TimeDiscretizationCompact( params.attrs['num_TimeSteps'], params.attrs['delta_t'] )

    return all, disc 



# =============== 12) Import CspinDMFT Data =============== (refactored code No. 2) IDENTICAL TO NL-SPINDMFT EXCEPT FOR FOLDERNAME
# everything is saved in a single hdf5 file
def ImportData_CspinDMFT( spin_model, config, other_physical_data = None, project = None, extension = None, warn_ETV = True, postfix = None ):
    # determine the folder:
    root_folder = "CspinDMFT/Data/"
    foldername = Orientation(postfix) + root_folder
    if project != None:
        foldername += project + "/"

    # interpret input:
    if isinstance(spin_model,str):
        spin_model = "spinmodel=" + spin_model
    else: # spin_model is a list of two elements
        spin_model = "spinmodel=" + spin_model[0] + "__mfspinmodel=" + spin_model[1]
    if other_physical_data != None:
        other_physical_data = "__" + other_physical_data
    else:
        other_physical_data = ""
    if extension != None:
        extension = "_" + extension
    else:
        extension = ""

    # create file str:
    filename = foldername + spin_model + "__config=" + config + other_physical_data + extension

    # check whether the parameters file exists and if required add suffix ETV (eigenvalue threshold violated):
    ETV = False
    if not exists( filename + ".hdf5" ): # is the negative eigenvalues threshold is violated in the data (then the foldername obtains an extension)
        filename = filename + "_ETV"
        ETV = True
        if not exists( filename + ".hdf5" ):
            raise Exception( "Could not find the file : %s" %(filename + ".hdf5") )
    filename += ".hdf5"
    all = h5.File( filename, 'r' )

    # warn if warner is active:
    if ETV and warn_ETV: 
        warnings.warn( "Eigenvalue threshold violated in data set of file \'%s\'" %filename )

    # discretization
    params =  all['parameters']
    disc = TimeDiscretizationCompact( params.attrs['num_TimeSteps'], params.attrs['delta_t'] )

    return all, disc 

def get_all_Gs( all ):
    return [[4*gab for gab in all['results']['correlation'][key]] for key in all['results']['correlation'].keys()]



# =============== 13) Import Classical spinDMFT Data ===============
# everything is saved in a single hdf5 file
def ImportData_Classical_spinDMFT( spin_model, physical_data = "", project = "", selfcons = True, extension = "" ):

    # process the inserted data:
    root_folder = "Classical_spinDMFT/Data/"
    if physical_data != "":
        physical_data = "__" + physical_data
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation() + root_folder
    if not selfcons:
        foldername += "noselfcons/"
    if project != "":
        foldername += project + "/"

    # determine file and return data:
    filename = foldername + "spinmodel=" + spin_model + physical_data + extension + ".hdf5"
    all = h5.File( filename, 'r' )

    # discretization
    disc = TimeDiscretizationCompact( all['parameters'].attrs['num_TimeSteps'], all['parameters'].attrs['delta_t'] )

    return all, disc


# =============== 14) Import p-spinDMFT Data ===============
# everything is saved in a single hdf5 file
def ImportData_pspinDMFT( spin_model, config, num_Sites, physical_data = "", project = "", selfcons = True, extension = "" ):

    # process the inserted data:
    root_folder = "p-spinDMFT/Data/"
    if physical_data != "":
        physical_data = "__" + physical_data
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = Orientation() + root_folder
    if not selfcons:
        foldername += "noselfcons/"
    if project != "":
        foldername += project + "/"

    # determine file and return data:
    filename = foldername + "spinmodel=" + spin_model + "__config=" + config + "__sites=" + str(num_Sites) + physical_data + extension + ".hdf5"
    all = h5.File( filename, 'r' )

    # discretization
    disc = TimeDiscretizationCompact( all['parameters'].attrs['num_TimeSteps'], all['parameters'].attrs['delta_t'] )

    return all, disc 

def get_g_pspinDMFT( all ):
    return [[gab for gab in all['results/correlation'][str(i)+"-"+str(i)]] for i in range(all['parameters'].attrs['num_Sites'])]
