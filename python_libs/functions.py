import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from io import StringIO
import os

# set labels:
def SetLabel(variable, unit, dim, style):
    if style=="":
        labelstr = r"$" + variable + r"\left[" + unit + r"\right]$"
    elif style=="PRB":
        labelstr = r"$" + variable + r"\,\left(\mathrm{units}\,\mathrm{of}\," + unit + r"\right)$"
    else: 
        raise Exception("unknown style")

    if dim=="x":
        plt.xlabel(r'%s'%labelstr)
    elif dim=="y":
        plt.ylabel(r'%s'%labelstr)
    else:
        raise Exception("unknown dimension")


# get labelstr:
def GetLabelStr(variable, unit, style):
    if style=="":
        labelstr = r"$" + variable + r"\left[" + unit + r"\right]$"
    elif style=="PRB":
        labelstr = r"$" + variable + r"\,\left(\mathrm{units}\,\mathrm{of}\," + unit + r"\right)$"
    else: 
        raise Exception("unknown style")

    return labelstr


# build a folder tree:
def BuildFolderTree( tree, path = "" ):
    if np.shape(tree) == ():
        total_path = str(tree)
        if path != "":
            if path[-1] != "/":
                total_path = path + "/" + total_path
            else:
                total_path = path + total_path
        try:
            os.mkdir( total_path )
        except:
            pass
    else:
        total_path = tree[0]
        if path != "":
            if path[-1] != "/":
                total_path = path + "/" + total_path
            else:
                total_path = path + total_path
        try:
            os.mkdir( total_path )
        except:
            pass
        for f in range( 1, len(tree) ):
            total_path += "/" + tree[f]
            try:
                os.mkdir( total_path )
            except:
                pass
    return total_path + "/"
            

# Curve Configuration:
class CURVECONFIG():
    # initial "private" variables:
    _color_counter = 0
    _marker_counter = 0
    _markevery_counter = 0
    _markerlimit = 8

    def __init__(self, numplots, cmap_name="inferno", cmap=0, color_steps="total", markermap_name="1", markermap=0, markeveryarray=np.array((40,44,48,52,56,60,64,68))): # number of plots (max. markerlimit), name of cmap, range of colors, name of markermap, markeveryarray
     
        # get numplots:
        self.numplots = numplots

        # get cmap:
        if cmap==0: # no cmap inserted => interpreting cmap_name
            self.cmap = plt.get_cmap('%s'%cmap_name)
        else:
            self.cmap = cmap
        
        # get markermap:
        if np.isscalar(markermap): # no markermap inserted => interpreting markermap_name
            if markermap_name=="1":
                self.markermap = np.array(('x','s','^','v','o','D','2','1','p','P','$u$'))
            else:
                raise Exception("unknown value for markermap_name")
        else:
            self.markermap=markermap

        # get markeveryarray:
        self.markeveryarray = markeveryarray

        # range of colors within cmap: || dc: step width, counter_plus: shift
        if color_steps=="left":
            self._dc = 1/numplots
            self._counter_plus = 0 
        elif color_steps=="right":
            self._dc = 1/numplots
            self._counter_plus = 1
        elif color_steps=="centered":
            self._dc = 1/numplots
            self._counter_plus = 0.5
        elif color_steps=="total":
            self._dc = 1/(numplots-1)
            self._counter_plus = 0
        else:
            raise Exception("unknown value for color_steps")

    def GetColor(self, increase=True):
        if self._color_counter==self.numplots: # maximum exceeded
            self._color_counter=0
            print("WARNING in class CURVECONFIG: Total color set consumed. Coloring restarts.")
        color = self.cmap( (self._color_counter + self._counter_plus) * self._dc )
        if increase:
            self._color_counter+=1
        return color 

    def GetMarker(self, increase=True):
        if self._marker_counter>=self.numplots:
            self._marker_counter=0
            print("WARNING in class CURVECONFIG: Total marker set consumed. Marking restarts.")
            
        if self._marker_counter>=len(self.markermap):
            self._marker_counter=0
            print("WARNING in class CURVECONFIG: Markerlimit exceeded. Marking restarts.")

        marker = self.markermap[self._marker_counter]
        if increase:
            self._marker_counter+=1
        return marker

    def GetMarkEvery(self, increase=True):
        if self._markevery_counter>=self.numplots:
            self._markevery_counter=0
            print("WARNING in class CURVECONFIG: Total markevery set consumed. Markevery'ing restarts.")
            
        if self._markevery_counter>=len(self.markeveryarray):
            self._markevery_counter=0
            print("WARNING in class CURVECONFIG: Markeverylimit exceeded. Markevery'ing restarts.")

        markevery = self.markeveryarray[self._markevery_counter]
        if increase:
            self._markevery_counter+=1
        return markevery


