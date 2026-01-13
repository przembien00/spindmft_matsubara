import numpy as np 
import matplotlib.pyplot as plt 
import import_data as id

spinmodels = ["DRF",["DRF","Ising"],"DRF"]
configs = ["Test_1","Test_2","Test_3"]
site_names_list = [['1-1','1-2'],['1-1','1-2'],['1-1','1-2']]
dir_names_list = [["xx","zz"],["xx","zz"],["xx","zz"]]
dir_indices_list = [[0,1],[0,1],[0,1]]

for spinmodel,config,site_names,dir_names,dir_indices in zip(spinmodels,configs,site_names_list,dir_names_list,dir_indices_list):
    # import Reference and Test data
    all_ref, disc_ref = id.ImportData_nlspinDMFT(spinmodel,config,project="Tests",extension="Reference_Data")
    G_ref_list = [[4*gab for gab in all_ref['results']['correlation'][pc]] for pc in site_names]
    all, disc = id.ImportData_nlspinDMFT(spinmodel,config,project="Tests",extension="Test_Data")
    G_list = [[4*gab for gab in all['results']['correlation'][pc]] for pc in site_names]

    # plot the data
    counter = 0
    colors = ["orange","blue","lime","black","red","cyan","magenta","yellow","purple"]
    for G_ref, G in zip(G_ref_list,G_list): # loop over ij
        for dir_name,dir_index in zip(dir_names,dir_indices): # loop over alpha beta
            label1, label2 = None, None
            if counter == 0:
                label1 = "Reference for " + "|".join(dir_names) + " & " + "|".join(site_names)
                label2 = "Test for " + "|".join(dir_names) + " & " + "|".join(site_names)
            plt.plot(disc_ref.t, G_ref[dir_index], color=colors[counter], ls="--", label=label1)
            plt.plot(disc.t, G[dir_index], color=colors[counter], ls="-", label=label2)
            counter += 1

    # finalize
    plt.title(config)
    plt.legend(loc='best')
    plt.show()

    '''
    fname = "nlspinDMFT_"+config+".pdf"
    plt.savefig(fname)
    '''