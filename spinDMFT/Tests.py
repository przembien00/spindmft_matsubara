import numpy as np 
import matplotlib.pyplot as plt
import sys
sys.path.append("../python_libs") 
import import_data as id

Tests = [str(i+1) for i in range(0,4)]
spinmodels = ["ISO","DRF","DRF","DRF"]
phys_datas = ["","","h=z_h_abs=1",""]
plot_correlations = [[0],[0,1],[0,1,2,3],[0,1,8]]
correlation_names = [["xx"],["xx","zz"],["xx","xy","yx","zz"],["xx","xy","zz"]]

for i,spinmodel,physdat,pcs,cnames in zip(Tests,spinmodels,phys_datas,plot_correlations,correlation_names):
    # import Reference and Test data
    all_ref, disc_ref = id.ImportData_spinDMFT(spinmodel,physical_data=physdat,project="Tests",extension="Reference_"+i)
    G_ref = [4*gab for gab in all_ref['results']['correlation']]
    all, disc = id.ImportData_spinDMFT(spinmodel,physical_data=physdat,project="Tests",extension="Test_"+i)
    G = [4*gab for gab in all['results']['correlation']]

    # plot the data
    colors = ["orange","blue","lime","black","red","cyan","magenta","yellow","purple"]
    for pc,cname,color in zip(pcs,cnames,colors):
        plt.plot(disc_ref.t, G_ref[pc], color=color, ls="--", label="Reference for "+cname)
        plt.plot(disc.t, G[pc], color=color, ls="-", label="Test for "+cname)

    # finalize
    plt.title("Test "+i)
    plt.legend(loc='best')
    plt.show()

    '''
    fname = "spinDMFT_Test_"+i+".pdf"
    plt.savefig(fname)
    '''