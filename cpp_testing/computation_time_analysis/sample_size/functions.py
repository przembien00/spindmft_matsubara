import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
import os.path

# ---- this function reads out the desired data ----
# root_folder describes the root directory of the data: "algorithm", "expval", "pulses"
# model_frame: "isotropic_lab", "isotropic_rot", "dipole_lab", "dipole_rot", "XXchain_lab"
# symmetry_type: "A", "B", "C", "D"
# phys_data: string with the data for C, B, timedep, ... including delimiter ','
# extension: any name extensions of the data folder
# --------------------------------------------------
# The matrix returns are transposed, i.e., t1 denotes the column and t2 denotes the row
def read_data(root_folder, model_frame, symmetry_type, phys_data, extension):
    # orientation in the folder tree:
    tree_climb = ""
    i=0 # preventing an infinite loop
    while not os.path.isfile("./%s.orientation"%tree_climb) and i<10: # climb the folder tree until the root folder is reached
        tree_climb = tree_climb + "../"
        i+=1

    # process the inserted data:
    root_folder = np.asarray(root_folder.split(","))
    symmetry_type = np.asarray(symmetry_type.split(","))
    phys_data = np.asarray(phys_data.split(","))
    if extension != "":
        extension = "_" + extension

    # return data from the TTI folder:
    if root_folder[0]=="selfcons_TTI" or root_folder[0]=="expval_TTI":
        if root_folder[0]=="selfcons_TTI":
            filepre = tree_climb + "DMFT_selfcons_TTI/Data/" + model_frame + "/Type" + symmetry_type[0]
        elif root_folder[0]=="expval_TTI":
            filepre = tree_climb + "DMFT_expval_TTI/Data/" + root_folder[1] + "/" + model_frame + "/Type" + symmetry_type[0]
        if symmetry_type[0] == "A":
            filepre = filepre + "/C=" + phys_data[0] + extension
            gxx = np.genfromtxt("%s/results.txt"%filepre, unpack = 'True')
            A = np.array((gxx))
        elif symmetry_type[0] == "B":
            filepre = filepre + "/C=" + phys_data[0] + extension
            gxx,gzz = np.genfromtxt("%s/results.txt"%filepre, unpack = 'True')
            A = np.array((gxx,gzz))
        elif symmetry_type[0] == "C":
            filepre = filepre + "/C=" + phys_data[0] + "_B=" + phys_data[1] + extension
            gxx,gxy,gzz = np.genfromtxt("%s/results.txt"%filepre, unpack = 'True')
            A = np.array((gxx,gxy,gzz))
        elif symmetry_type[0] == "D":
            filepre = filepre + "/C=" + phys_data[0] + "_B=" + phys_data[1] + "_tau=" + phys_data[2] + "_varphi=" + phys_data[3] + extension
            gxx,gxy,gxz,gyx,gyy,gyz,gzx,gzy,gzz = np.genfromtxt("%s/results.txt"%filepre, unpack = 'True')
            A = np.array((gxx,gxy,gxz,gyx,gyy,gyz,gzx,gzy,gzz))
        steps, Delta_t, M = np.genfromtxt("%s/params.txt"%filepre, unpack = 'True')
        t = np.linspace(0, steps*Delta_t, steps+1)
        return A, t

    # return data from the PUL folder:
    elif root_folder[0]=="selfcons_PUL":
        filepre = tree_climb + "DMFT_selfcons_PUL/Data/" + model_frame + "/Type" + symmetry_type[0]
        if symmetry_type[0] == "A":
            filepre = filepre + "/B=0.00" + extension
            gxx = np.loadtxt("%s/results_gxx.txt"%filepre)
            A = np.array((gxx))
        elif symmetry_type[0] == "B":
            filepre = filepre + "/B=0.00" + extension
            gxx = np.genfromtxt("%s/results_gxx.txt"%filepre, unpack = 'True')
            gzz = np.genfromtxt("%s/results_gzz.txt"%filepre, unpack = 'True')
            A = np.array((gxx,gzz))
        elif symmetry_type[0] == "C":
            filepre = filepre + "/tdp=" + phys_data[0] + "_B=" + phys_data[1] + "_param1=" + phys_data[2] + "_param2=" + phys_data[3] + extension
            gxx = np.genfromtxt("%s/results_gxx.txt"%filepre, unpack = 'True')
            gxy = np.genfromtxt("%s/results_gxy.txt"%filepre, unpack = 'True')
            gzz = np.genfromtxt("%s/results_gzz.txt"%filepre, unpack = 'True')
            A = np.array((gxx,gxy,gzz))
        elif symmetry_type[0] == "D":
            filepre = filepre + "/tdp=" + phys_data[0] + "_B=" + phys_data[1] + "_param1=" + phys_data[2] + "_param2=" + phys_data[3] + extension
            gxx = np.genfromtxt("%s/results_gxx.txt"%filepre, unpack = 'True')
            gxy = np.genfromtxt("%s/results_gxy.txt"%filepre, unpack = 'True')
            gxz = np.genfromtxt("%s/results_gxz.txt"%filepre, unpack = 'True')
            gyy = np.genfromtxt("%s/results_gyy.txt"%filepre, unpack = 'True')
            gyz = np.genfromtxt("%s/results_gyz.txt"%filepre, unpack = 'True')
            gzz = np.genfromtxt("%s/results_gzz.txt"%filepre, unpack = 'True')
            A = np.array((gxx,gxy,gxz,gyy,gyz,gzz))
        steps, Delta_t, M = np.genfromtxt("%s/params.txt"%filepre, unpack = 'True')
        t = np.linspace(0, steps*Delta_t, steps+1)
        return A, t
    elif root_folder[0]=="expval_PUL":
        filepre = tree_climb + "DMFT_expval_PUL/Data/" + root_folder[1] + "/" + model_frame + "/Type" + symmetry_type[0]
        src_data = ""
        l = len(phys_data)-1
        if symmetry_type[1] == "A":
            src_data = src_data + "_B=0.00"
        elif symmetry_type[1] == "B":
            src_data = src_data + "_B=0.00"
        elif symmetry_type[1] == "C":
            src_data = src_data +  "_tdp=" + phys_data[l-3] + "_B=" + phys_data[l-2] + "_param1=" + phys_data[l-1] + "_param2=" + phys_data[l]
        elif symmetry_type[1] == "D":
            src_data = src_data +  "_tdp=" + phys_data[l-3] + "_B=" + phys_data[l-2] + "_param1=" + phys_data[l-1] + "_param2=" + phys_data[l]
        if symmetry_type[0] == "A":
            filepre = filepre + "/B=0.00" + extension
            gxx = np.loadtxt("%s/results_gxx.txt"%filepre)
            A = np.array((gxx))
        elif symmetry_type[0] == "B":
            filepre = filepre + "/B=0.00" + extension
            gxx = np.genfromtxt("%s/results_gxx.txt"%filepre, unpack = 'True')
            gzz = np.genfromtxt("%s/results_gzz.txt"%filepre, unpack = 'True')
            A = np.array((gxx,gzz))
        elif symmetry_type[0] == "C":
            filepre = filepre + "/pulse=" + phys_data[0] + "_pB=" + phys_data[1] + "_pparam1=" + phys_data[2] + "_pparam2=" + phys_data[3] + src_data + extension
            gxx = np.genfromtxt("%s/results_gxx.txt"%filepre, unpack = 'True')
            gxy = np.genfromtxt("%s/results_gxy.txt"%filepre, unpack = 'True')
            gzz = np.genfromtxt("%s/results_gzz.txt"%filepre, unpack = 'True')
            A = np.array((gxx,gxy,gzz))
        elif symmetry_type[0] == "D":
            filepre = filepre + "/pulse=" + phys_data[0] + "_pB=" + phys_data[1] + "_pparam1=" + phys_data[2] + "_pparam2=" + phys_data[3] + src_data + extension
            gxx = np.genfromtxt("%s/results_gxx.txt"%filepre, unpack = 'True')
            gxy = np.genfromtxt("%s/results_gxy.txt"%filepre, unpack = 'True')
            gxz = np.genfromtxt("%s/results_gxz.txt"%filepre, unpack = 'True')
            gyy = np.genfromtxt("%s/results_gyy.txt"%filepre, unpack = 'True')
            gyz = np.genfromtxt("%s/results_gyz.txt"%filepre, unpack = 'True')
            gzz = np.genfromtxt("%s/results_gzz.txt"%filepre, unpack = 'True')
            A = np.array((gxx,gxy,gxz,gyy,gyz,gzz))
        steps, Delta_t, M = np.genfromtxt("%s/params.txt"%filepre, unpack = 'True')
        t = np.linspace(0, steps*Delta_t, steps+1)
        return A, t

def param_file_to_strlist(root_folder, model_frame, symmetry_type, phys_data, extension):
    # orientation in the folder tree:
    tree_climb = ""
    i=0 # preventing an infinite loop
    while not os.path.isfile("./%s.orientation"%tree_climb) and i<10: # climb the folder tree until the root folder is reached
        tree_climb = tree_climb + "../"
        i+=1

    # process the inserted data:
    root_folder = np.asarray(root_folder.split(","))
    symmetry_type = np.asarray(symmetry_type.split(","))
    phys_data = np.asarray(phys_data.split(","))
    if extension != "":
        extension = "_" + extension

    filepre = tree_climb + "DMFT_selfcons_PUL/Data/" + model_frame + "/Type" + symmetry_type[0]
    if symmetry_type[0] == "A":
        filepre = filepre + "/B=0.00" + extension
    elif symmetry_type[0] == "B":
        filepre = filepre + "/B=0.00" + extension
    elif symmetry_type[0] == "C":
        filepre = filepre + "/tdp=" + phys_data[0] + "_B=" + phys_data[1] + "_param1=" + phys_data[2] + "_param2=" + phys_data[3] + extension
    elif symmetry_type[0] == "D":
        filepre = filepre + "/tdp=" + phys_data[0] + "_B=" + phys_data[1] + "_param1=" + phys_data[2] + "_param2=" + phys_data[3] + extension
    with open ("%s/params.txt"%filepre, "r") as myfile:
        data = myfile.readlines()
    return data

# returns the position of correlation corr in the return array A:
def plot_number(corr):
    if corr == "xx":
        return 0,0,0,0
    elif corr == "xy":
        return -1,-1,1,1
    elif corr == "xz":
        return -1,-1,-1,2
    elif corr == "yy":
        return 0,0,0,3
    elif corr == "yz":
        return -1,-1,-1,4
    elif corr == "zz":
        return 0,1,2,5
