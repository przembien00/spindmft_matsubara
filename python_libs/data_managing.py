import numpy as np 
import os
import inspect
import h5py as h5


DATAMANAGER_description = "manages datasets created in python\n" 
class DATAMANAGER():
    # initialization:
    def __init__(self,save_as='hdf5'):
        self.names = []
        self.datasets = []
        self.descriptions = []
        self.infos = []
        self.save_as = save_as

    # list all available methods:
    def help(self):
        print( "========== class DATAMANAGER ==========\n")
        print( "=== description ===" )
        print( DATAMANAGER_description )
        print( "=== available methods ===")
        for i in inspect.getmembers(self):
            if not i[0].startswith('_'): # remove private and protected stuff
                if inspect.ismethod(i[1]):
                    print(i[0])

    # extract and locally store a dataset:
    def extract_dataset(self,name,data,description,**kwargs):
        self.names.append(name)
        self.datasets.append(data)
        self.descriptions.append(description)
        self.infos.append(kwargs)

    # store the data:
    def store(self,filename):
        if self.save_as == 'hdf5':
            with h5.File(filename+'.'+self.save_as,'w') as file:
                for name,dataset,description,info in zip(self.names,self.datasets,self.descriptions,self.infos):
                    dset = file.create_dataset(name,data=dataset)
                    dset.attrs['description'] = description
                    for attribute in info:
                        dset.attrs[attribute[0]] = attribute[1]
        else:
            raise RuntimeError("storing format %s undefined"%self.save_as)