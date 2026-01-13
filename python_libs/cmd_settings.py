import warnings as warn
import numpy as np 
import sys 

CRED = '\033[91m'
CEND = '\033[0m'

# interprets the command line arguments and returns them in the appropriate data type
def cmdline_parameters( datatypes ):
    parameters = sys.argv[1:]
    datatypes = np.asarray( datatypes.split(',') )
    if len(parameters) != len(datatypes): # raise error
        raise RuntimeError( "number of datatypes doesn't match the number of cmd parameters" )
    else: # return cmd arguments
        for index in range(0, len(parameters)):
            parameters[index] = transform_to_type( parameters[index], datatypes[index] )
        return parameters

# same as above, but with default values, if number oc cmd parameters does not match the required number
def cmdline_parameters_with_defaults( datatypes, defaults ):
    parameters = sys.argv # first is filename
    datatypes = np.asarray( datatypes.split(',') )
    defaults = np.asarray( defaults.split(',') )
    if len(parameters)-1 != len(datatypes): # return the default values
        print(CRED, "Using default parameters in this run", CEND)
        for index in range(0, len(defaults)):
            defaults[index] = transform_to_type( defaults[index], datatypes[index] )
        return defaults
    else: # return cmd arguments
        parameters = parameters[1:] # ignore filename
        for index in range(0, len(parameters)):
            parameters[index] = transform_to_type( parameters[index], datatypes[index] )
        return parameters

# transforms a parameter into a desired type
def transform_to_type( p, dtype ):
    if dtype=="int":
        return int(p)
    elif dtype=="flt":
        return float(p)
    elif dtype=="str":
        return str(p)
    else:
        raise RuntimeError( "dtype not implemented in transform_to_type" )
