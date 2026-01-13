import numpy as np 
from scipy.signal import find_peaks


# ================================================================
# ======================== PROCESS CURVES ========================
# ================================================================


# function : computes the envelopes of some oscillating data by finding minima and maxima
def ComputeEnvelope( x, signal ):
    numsteps = len(signal)
    upper_env = [] 
    upper_env_x = []
    lower_env = []
    lower_env_x = []

    # find the minima and maxima and add them to the envelopes
    # using FindMaxima or FindMinima would be inefficient here
    for i in range( 1, numsteps-1 ):
        if signal[i] > signal[i-1]: # potential maxima
            if signal[i] > signal[i+1]: # maxima
                upper_env = upper_env + [signal[i]]
                upper_env_x = upper_env_x + [x[i]]
        else: # potential minima
            if signal[i] < signal[i+1]: # minima
                lower_env = lower_env + [signal[i]]
                lower_env_x = lower_env_x + [x[i]]
    
    return np.array(lower_env_x), np.array(lower_env), np.array(upper_env_x), np.array(upper_env)


# function : computes the Nth envelopes of some oscillating data by finding minima and maxima 
# this can be used for weak discretizations that do not fully resolve the oscillations
def ComputeNthEnvelope( x, signal, N ):
    lower_env_x, lower_env, upper_env_x, upper_env = ComputeEnvelope( x, signal )

    for n in range(1,N): # compute the envelope of the envelope N-1 times
        lower_env_x, lower_env, waste_1, waste_2 = ComputeEnvelope( lower_env_x, lower_env )
        waste_1, waste_2, upper_env_x, upper_env = ComputeEnvelope( upper_env_x, upper_env )
    
    return np.array(lower_env_x), np.array(lower_env), np.array(upper_env_x), np.array(upper_env)


# function : computes the avergage of the upper and lower Nth envelopes
# assumes x is sorted
def ComputeNthEnvelopeAverage( x, signal, N ):
    # compute the Nth envelope
    lower_env_x, lower_env, upper_env_x, upper_env = ComputeNthEnvelope( x, signal, N )

    # choose a new equidistant discretization
    dx = ( lower_env_x[-1] + lower_env_x[0] ) / len(lower_env_x) # step width
    x_min = np.max( [lower_env_x[0],upper_env_x[0]] ) # lower bound
    x_max = np.min( [lower_env_x[-1],upper_env_x[-1]] ) - dx/2 - dx/4 # upper bound
    x_new = [x_min + dx/2]
    while x_new[-1] < x_max:
        x_new = x_new + [x_new[-1] + dx/2]

    # extrapolate the envelopes to this discretization
    waste, new_lower_env = ExtrapolateLinearly( lower_env_x, lower_env, x_new )
    waste, new_upper_env = ExtrapolateLinearly( upper_env_x, upper_env, x_new )

    # average the envelopes 
    envelope_average = (new_lower_env + new_upper_env)/2
    
    return x_new, envelope_average


# function : extrapolates a signal linearly to a certain domain point x_m
# assumes x is sorted
def ExtrapolateDataPointLinearly( x, signal, x_m ):
    # find the closest values to x_m in x
    index_closest = np.argmin(np.abs(x - x_m))

    # determine the indices of the surrounding points if possible
    if x[index_closest] > x_m: # above
        if index_closest != 0: # closest value is not the first value
            index_up  = index_closest
            index_low = index_closest-1
        else: # closest value is above x_m and first value at the same time -> no extrapolation
            return None, None
    else: # below
        if index_closest != len(x)-1: # closest value is not the last value
            index_up  = index_closest+1
            index_low = index_closest
        else: # closest value is below x_m and last value at the same time -> no extrapolation
            return None, None
    
    # linear fit         
    x_1 = x[index_low]
    x_2 = x[index_up]
    y_1 = signal[index_low]
    y_2 = signal[index_up]
    m = (y_2-y_1) / (x_2-x_1)
    b = (y_1*x_2-y_2*x_1) / (x_2-x_1)
    y_m = m*x_m + b

    return x_m, y_m


# function : extrapolates a signal linearly to a new discretization
# automatically reduces x_new if it goes beyond the limits
# assumes x and x_new are sorted
# in case x and x_new start at the same value the first point is treated separately
def ExtrapolateLinearly( x, signal, x_new ): 
    new_signal = []
    new_x = []
    if x[0] == x_new[0]: # then treat the first point separately
        new_signal = new_signal + [signal[0]]
        new_x = new_x + [x[0]]
        start = 1
    else:
        start = 0

    for i in range(start, len(x_new)):
        x_m, y_m = ExtrapolateDataPointLinearly(x, signal, x_new[i])
        if not x_m or not y_m: # no extrapolation possible
            continue

        # add to data 
        new_x = new_x + [x_m]
        new_signal = new_signal + [y_m]

    return np.array(new_x), np.array(new_signal)

# Eliminate all but the nth data points of a signal
def ReduceDiscretization( data, n ):
    reduced_data = []
    for i in range(0, len(data)):
        if i % n == 0:
            reduced_data = reduced_data + [data[i]]
    return np.array( reduced_data )

# adds two signals together
# assumes both signals have the same discretization but different amount of data points 
# the excess data points are cutted
def AddAndCut( signal1, signal2 ):
    new_signal = []
    for s1,s2 in zip(signal1,signal2): # indices behind the last joint one are ignored
        new_signal.append(s1+s2)
    return np.array(new_signal)

# rescale a signal several times with different rescaling factors
# return list of rescaled signals
def RescaleToList( x, signal, J0, Js ):
    resc_signals = []
    for J in Js:
        new_x = x * J0/J
        trash, resc_signal = ExtrapolateLinearly(new_x,signal,x)
        resc_signals.append(resc_signal)
    return resc_signals


# accumulates a list of signals (arrays)
def AccumulateAndCut( signals, weights=() ):
    max_signal_length = np.min([len(signal) for signal in signals])
    acc_signal = []
    if np.shape(weights) == (0,):
        for i in range(0,max_signal_length):
            signals_at_same_time = [s[i] for s in signals]
            acc_signal.append( np.sum(signals_at_same_time) )
    elif np.shape(weights) == ():
        weight = weights
        for i in range(0,max_signal_length):
            signals_at_same_time = [weight*s[i] for s in signals]
            acc_signal.append( np.sum(signals_at_same_time) )
    else:
        for i in range(0,max_signal_length):
            signals_at_same_time = [weight*s[i] for s,weight in zip(signals,weights)]
            acc_signal.append( np.sum(signals_at_same_time) )
    return acc_signal

# accumulates a list of signals (arrays)
def AverageAndCut( signals ):
    num_signals = len( signals )
    return AccumulateAndCut( signals, weights=1/num_signals )

# function : computes the square difference of two signals and normalizes it
# assumes both signals on the same discretization
def ComputeSquareAveragedError( signal1, signal2 ):
    error = 0
    for s1,s2 in zip(signal1,signal2):
        error += (s1-s2)**2
    return error/len(signal1)

# determine optimal rescaling for a dataset (d2,s2) with respect to another dataset (d1,s1)
# the best match is searched considering the squared average error
def DetermineRescale( d1, s1, d2, s2, rescale_region, grid_points = 1000 ):
    rescale = np.linspace(rescale_region[0], rescale_region[1], grid_points)
    error = np.zeros(len(rescale))
    for index,r in enumerate(rescale):
        d2_new = r * d2 
        d_new, s2_new = ExtrapolateLinearly( d2_new, s2, d1 ) # extrapolate s2 to d1
        first = np.argwhere( d1 == d_new[0] )[0,0]
        last = np.argwhere( d1 == d_new[-1] )[0,0]
        s1_new = s1[first:last+1]
        error[index] = ComputeSquareAveragedError( s2_new, s1_new )
    return rescale[np.argmin(error)]


# function : computes the time discretization error over time (fine signal - coarse signal)
# assumes both signals are defined on the same time interval
# assumes a power of 2 as the ratio between the number of time steps 
def DiscretizationError_pow2( signal_fine, signal_coarse ):
    ratio = int((len(signal_fine)-1)/(len(signal_coarse)-1)) # ratio of the number of time steps

    # discretization error over time
    dq_over_time = np.zeros(len(signal_coarse))
    for i in range(0, len(dq_over_time)):
        dq_over_time[i] = signal_fine[i*ratio] - signal_coarse[i]

    # return time average
    return dq_over_time

# function : computes the deviation between two signals
# assumes both signals are defined on the same time discretization
def Deviation( signal_1, signal_2 ):
    # error over time
    dq_over_time = []
    for s1, s2 in zip(signal_1, signal_2):
        dq_over_time.append(s1 - s2)

    # return time average
    return dq_over_time


# function : computes the time-averaged time discretization error 
# assumes both signals are defined on the same time interval
# assumes a power of 2 as the ratio between the number of time steps 
def AverageDiscretizationError_pow2( signal_fine, signal_coarse ):
    dq_over_time = DiscretizationError_pow2( signal_fine, signal_coarse )
    dq = np.mean(dq_over_time)
    dq_std = np.std(dq_over_time)
    return dq, dq_std


# function : finds the maxima of a signal and returns the indices
def FindMaxima( signal ):
    max_positions = []
    for i in range( 1, len(signal)-1 ):
        if signal[i] > signal[i-1]: # potential maxima
            if signal[i] > signal[i+1]: # local maxima at i 
                max_positions = max_positions + [i]
    return max_positions


# function : finds the minima of a signal and returns the indices
def FindMinima( signal ):
    return FindMaxima( -signal )



# LISTS
# create empty 1D list
def CreateEmpty2DList( size ):
    list = []
    for i in range(0,size):
        list.append([])
    return list

# create empty 2D list
def CreateEmpty3DList( size, subsize ):
    list = []
    for i in range(0,size):
        list.append(CreateEmpty2DList(subsize))
    return list

# return the max indices of a 2D list as tuple
def max_indices_of_2D_list(llist):
    local_max_indices = []
    local_max_vals = []
    for sublistindex, sublist in enumerate(llist):
        if sublist:
            max_val = np.max(sublist)
            indices = np.where(sublist==max_val)[0]
            for index in indices:
                local_max_indices.append([sublistindex,index])
                local_max_vals.append(max_val)
    if local_max_indices:
        max = np.max(local_max_vals)
        indices = np.where(local_max_vals==max)[0]
        max_indices = []
        for index in indices:
            max_indices.append(local_max_indices[index])
        return max_indices
    else:
        return []



# STRINGS AND PRINTING
# nice print function for curve fit results
def print_pparams( popt, pcov, names = "default", round = lambda x : np.around(x,8) ):
    if names == "default":
        names = [str(i) for i in range(1,len(popt)+1)]
    for po, pc, name, in zip(popt,np.diag(pcov),names):
        print(name, ":", round(po), round(np.sqrt(pc)))
    
# thanks chatgpt
def remove_trailing_zeros(number_string):
    # Convert the string to a float to handle cases like "003.00"
    number = float(number_string)
    # Convert the float back to a string without trailing zeros
    result = str(number).rstrip('0').rstrip('.') if '.' in str(number) else str(int(number))
    return result