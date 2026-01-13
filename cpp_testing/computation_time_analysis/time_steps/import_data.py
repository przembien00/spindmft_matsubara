import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
import os.path

# =============== reading out the data: ===============
# root_folder   : "selfcons", "expval"
# model_frame   : "isotropic_lab", "dipole_lab",...
# symmetry_type : "A", "B", "C", "D"
# phys_data     : string with the physical data including delimiter ','
# extension     : any name extensions of the data folder
# --------------------------------------------------
# The matrix returns are transposed, i.e., t1 denotes the column and t2 denotes the row

# PUL data
def read_PUL_data(root_folder, model_frame, symmetry_type, phys_data, extension):

    # process the inserted data:
    root_folder = np.asarray(root_folder.split(","))
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = orientation()
    if root_folder[0]=="selfcons":
        foldername += "DMFT_selfcons_PUL/Data/" + model_frame + "/Type" + symmetry_type + "/" + phys_data + extension
    elif root_folder[0]=="expval":
        foldername += "DMFT_expval_PUL/Data/" + root_folder[1] + "/" + model_frame + "/Type" + symmetry_type + "/" + phys_data + extension

    # return the data:
    data, params = read_PUL(foldername,symmetry_type)
    return data, return_time(foldername)

def read_PUL_params(root_folder, model_frame, symmetry_type, phys_data, extension):

    # process the inserted data:
    root_folder = np.asarray(root_folder.split(","))
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = orientation()
    if root_folder[0]=="selfcons":
        foldername += "DMFT_selfcons_PUL/Data/" + model_frame + "/Type" + symmetry_type + "/" + phys_data + extension
    elif root_folder[0]=="expval":
        foldername += "DMFT_expval_PUL/Data/" + root_folder[1] + "/" + model_frame + "/Type" + symmetry_type + "/" + phys_data + extension

    # return the data:
    data, params = read_PUL(foldername,symmetry_type)
    return params


# TTI data
def read_TTI_data(root_folder, model_frame, symmetry_type, phys_data, extension):

    # process the inserted data:
    root_folder = np.asarray(root_folder.split(","))
    if extension != "":
        extension = "_" + extension

    # determine the folder:
    foldername = orientation()
    if root_folder[0]=="selfcons":
        foldername += "DMFT_selfcons_TTI/Data/" + model_frame + "/Type" + symmetry_type + "/" + phys_data + extension
    elif root_folder[0]=="expval":
        foldername += "DMFT_expval_TTI/Data/" + root_folder[1] + "/" + model_frame + "/Type" + symmetry_type + "/" + phys_data + extension

    # return the data:
    data, params = read_PUL(foldername,symmetry_type)
    return data, return_time(foldername)

def make_array(phys_data):
    f = 0
    t = 0
    i = 0
    while t==0 and i<len(phys_data):
        if phys_data[i] == "[":
            f = i
        if phys_data[i] == "]":
            t = i
        i = i+1
    if f!=0 and t!=0:
        params = np.asarray(phys_data[f+1:t].split(','))
        strings = ["" for x in range(len(params))]
        for l in range(0,len(params)):
            strings[l] = phys_data[0:f] + params[l] + phys_data[t+1:]
        phys_data = strings
    else:
        return phys_data
    return phys_data

# returns the position of correlation corr in the return array A:
def plot_number(corr):
    if corr == "xx":
        return 0,0,0,0
    elif corr == "xy":
        return float('NaN'),float('NaN'),1,1
    elif corr == "xz":
        return float('NaN'),float('NaN'),float('NaN'),2
    elif corr == "yy":
        return 0,0,0,3
    elif corr == "yz":
        return float('NaN'),float('NaN'),float('NaN'),4
    elif corr == "zz":
        return 0,1,2,5

def orientation(): # orientation in the folder tree
    tree_climb = ""
    i=0 # preventing an infinite loop
    while not os.path.isfile("./%s.orientation"%tree_climb) and i<10: # climb the folder tree until the root folder is reached
        tree_climb = tree_climb + "../"
        i+=1
    return tree_climb

def return_time(foldername):
    steps, Delta_t, M = np.genfromtxt("%s/params.txt"%foldername, unpack = 'True')
    return np.linspace(0, steps*Delta_t, steps+1)

def read_PUL(foldername, symmetry_type):
    with open ("%s/params.txt"%foldername, "r") as myfile:
        params = myfile.readlines()
    if symmetry_type == "A":
        gxx = np.genfromtxt("%s/results_gxx.txt"%foldername, unpack='True')
        return np.array((gxx)), params
    elif symmetry_type == "B":
        gxx = np.genfromtxt("%s/results_gxx.txt"%foldername, unpack='True')
        gzz = np.genfromtxt("%s/results_gzz.txt"%foldername, unpack='True')
        return np.array((gxx,gzz)), params
    elif symmetry_type == "C":
        gxx = np.genfromtxt("%s/results_gxx.txt"%foldername, unpack = 'True')
        gxy = np.genfromtxt("%s/results_gxy.txt"%foldername, unpack = 'True')
        gzz = np.genfromtxt("%s/results_gzz.txt"%foldername, unpack = 'True')
        return np.array((gxx,gxy,gzz)), params
    elif symmetry_type == "D":
        gxx = np.genfromtxt("%s/results_gxx.txt"%foldername, unpack = 'True')
        gxy = np.genfromtxt("%s/results_gxy.txt"%foldername, unpack = 'True')
        gxz = np.genfromtxt("%s/results_gxz.txt"%foldername, unpack = 'True')
        gyy = np.genfromtxt("%s/results_gyy.txt"%foldername, unpack = 'True')
        gyz = np.genfromtxt("%s/results_gyz.txt"%foldername, unpack = 'True')
        gzz = np.genfromtxt("%s/results_gzz.txt"%foldername, unpack = 'True')
        return np.array((gxx,gxy,gxz,gyy,gyz,gzz)), params
