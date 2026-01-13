import matplotlib.pyplot as plt 
import numpy as np
import matplotlib.gridspec as gridspec

# MULTIPLOT WITH 2 COLUMNS AND N ROWS, X-AXIS SHARED
# dense plot with adjustable wspace and labels on the edges
# n: number of plots vertically, wspace_ratio: horizontal space between the plots, ratio_left: ratio of the left plot of the plotdwidth ("rightsquare" -> right plots are quadratic)
def dense_2column_plot( n = 2, wspace_ratio = 0.004, ratio_left = 0.5, ticks = "everywhere", height = 1.0 ): 
    
    # general stuff:
    h_0 = 1.89 * height # norm height of a plot
    w_0 = 5.67 / 2 # norm width of a plot
    margin = 0.71  # norm margin
    fwidth = 2*w_0 + 2*margin # figure width  || scales to 18 cm
    fheight = n*h_0 + margin  # figure height || scales to 4.80 cm + 1.8 cm (for height = 1)
    fheight_0 = h_0 + margin

    # edges of the plot boxes:
    t_e = 1 - 0.05 * fheight_0 / fheight # scaled top edge
    b_e = ( margin / fheight_0 - 0.05 ) * fheight_0 / fheight  # scaled bottom edge
    l_e = 0.1
    r_e = 0.9  

    # horizontal space between the plots:
    wspace = wspace_ratio * (r_e - l_e) * w_0 # scales to 0.23 cm * wspace_ratio

    # ratio of the left plot width from the whole plots width:
    if ratio_left == "rightsquare": # make the right plot squared
        ratio_left = ( 2*w_0 - h_0 ) / ( 2*w_0 )

    # generate figure and plotgrid:
    fig = plt.figure( figsize = (fwidth,fheight) )
    gs = gridspec.GridSpec( n, 2, wspace=wspace, hspace=0, top=t_e, bottom=b_e, left=l_e, right=r_e, width_ratios=(ratio_left,1-ratio_left) )

    # generate axis array:
    ax = np.empty((n,2), dtype=object)
    ax[n-1,0] = plt.subplot( gs[2*n-2] ) # first the bottom plot to share the axis 
    ax[n-1,1] = plt.subplot( gs[2*n-1] ) # first the bottom plot to share the axis
    for axis in range(0,n-1): # then the other axes
        ax[axis,0] = plt.subplot( gs[axis*2]  , sharex = ax[n-1,0] ) # share axis
        ax[axis,1] = plt.subplot( gs[axis*2+1], sharex = ax[n-1,1] ) # share axis
        ax[axis,0].tick_params( 'x', labelbottom = False ) # hide ticks 
        ax[axis,1].tick_params( 'x', labelbottom = False ) # hide ticks

    # set ticks:
    if ticks == "everywhere":
        for axis in range(0,n):
            ax[axis,1].yaxis.tick_right() # ticks to the right
            ax[axis,0].tick_params( right = True )
            ax[axis,1].tick_params( left = True )
            ax[axis,0].tick_params( right = True, which = 'minor' )
            ax[axis,1].tick_params( left = True, which = 'minor'  )
    elif ticks == "half":
        for axis in range(0,n):
            ax[axis,1].yaxis.tick_right() # ticks to the right
    else:
        raise Exception("tick setting %s unknown"%ticks)

    # set labels:
    for axis in range(0,n):
        ax[axis,1].yaxis.set_label_position("right") # labels to the right

    return ax


# MULTIPLOT WITH 2 COLUMNS AND N ROWS, NO AXIS SHARED
# dense plot with adjustable wspace and labels on the edges
# n: number of plots vertically, wspace_ratio: horizontal space between the plots, ratio_left: ratio of the left plot of the plotdwidth ("rightsquare" -> right plots are quadratic)
def dense_2column_plot_noshare( n = 2, wspace_ratio = 0.004, hspace_ratio = 0.2, ratio_left = 0.5, ticks = "everywhere", height = 1.0 ): 
    
    # general stuff:
    h_0 = 1.89 * height # norm height of a plot
    w_0 = 5.67 / 2 # norm width of a plot
    margin = 0.71  # norm margin
    fwidth = 2*w_0 + 2*margin # figure width  || scales to 18 cm
    fheight = n*h_0 + margin + (n-1) * margin * hspace_ratio # figure height || scales to 4.80 cm + 1.8 cm (for height = 1)
    fheight_0 = h_0 + margin

    # edges of the plot boxes:
    t_e = 1 - 0.05 * fheight_0 / fheight # scaled top edge
    b_e = ( margin / fheight_0 - 0.05 ) * fheight_0 / fheight  # scaled bottom edge
    l_e = 0.1
    r_e = 0.9  

    # horizontal space between the plots:
    wspace = wspace_ratio * (r_e - l_e) * w_0 # scales to 0.23 cm * wspace_ratio

    # vertical space between the plots:
    hspace = hspace_ratio * margin 

    # ratio of the left plot width from the whole plots width:
    if ratio_left == "rightsquare": # make the right plot squared
        ratio_left = ( 2*w_0 - h_0 ) / ( 2*w_0 )

    # generate figure and plotgrid:
    plt.figure( figsize = (fwidth,fheight) )
    gs = gridspec.GridSpec( n, 2, wspace=wspace, hspace=hspace, top=t_e, bottom=b_e, left=l_e, right=r_e, width_ratios=(ratio_left,1-ratio_left) )

    # generate axis array:
    ax = np.empty((n,2), dtype=object)
    for axis in range(0,n): # then the other axis
        ax[axis,0] = plt.subplot( gs[axis*2]   ) # share axis
        ax[axis,1] = plt.subplot( gs[axis*2+1] ) # share axis

    # set ticks:
    if ticks == "everywhere":
        for axis in range(0,n):
            ax[axis,1].yaxis.tick_right() # ticks to the right
            ax[axis,0].tick_params( right = True )
            ax[axis,1].tick_params( left = True )
            ax[axis,0].tick_params( right = True, which = 'minor' )
            ax[axis,1].tick_params( left = True, which = 'minor'  )
    elif ticks == "half":
        for axis in range(0,n):
            ax[axis,1].yaxis.tick_right() # ticks to the right
    else:
        raise Exception("tick setting %s unknown"%ticks)

    # set labels:
    for axis in range(0,n):
        ax[axis,1].yaxis.set_label_position("right") # labels to the right

    return ax



# MULTIPLOT WITH 2 COLUMNS AND N ROWS : Add subplotlabels
def dense_2column_plot_sublabels( n, ax, shift_x = [0,0], shift_y = [0,0] ): # number of subplot rows, axes, shift label in x/y dir (percentage)
    left_   = np.array(( '(a)', '(c)', '(e)', '(g)', '(i)', '(k)' ))
    right_  = np.array(( '(b)', '(d)', '(f)', '(h)', '(j)', '(l)' ))
    axis_   = range(0,n)
    for (left, right, axis) in zip(left_, right_, axis_):
        ax[axis,0].text(0.05 + shift_x[0], 0.95 + shift_y[0], left,  transform=ax[axis,0].transAxes, fontsize='medium', verticalalignment='top', fontfamily='serif' )
        ax[axis,1].text(0.05 + shift_x[1], 0.95 + shift_y[1], right, transform=ax[axis,1].transAxes, fontsize='medium', verticalalignment='top', fontfamily='serif' )



# MULTIPLOT WITH 2 COLUMNS AND N ROWS, X-AXIS SHARED
# dense plot with adjustable wspace and labels on the edges
# n: number of plots vertically, wspace_ratio: horizontal space between the plots, ratio_left: ratio of the left plot of the plotdwidth ("rightsquare" -> right plots are quadratic)
def dense_2column( f0size = (3, 2.1), rows = 2, wspace_ratio = 0.01, hspace_ratio = 0, ratio_left = 0.5, ticks = "everywhere", sharey = False ): 
    
    # general stuff:
    margin = 0.08*f0size[0] # norm horizontal and vertical margin
    fheight = rows*(f0size[1]-margin) + margin # 1 margin at the bottom
    fwidth = 1.8*f0size[0] # width of the plot (including 2 margins)
    wspace = wspace_ratio * fwidth # horizontal space between the plots
    hspace = hspace_ratio * fheight # horizontal space between the plots

    # ratio of the left plot width from the whole plots width:
    if ratio_left == "rightsquare": # make the right plot approximately squared
        ratio_right = (f0size[1]-margin)/(fwidth-2*margin)
        ratio_left = 1-ratio_right
    ratio_right = 1-ratio_left

    # generate figure and plotgrid:
    fig = plt.figure( figsize = (fwidth,fheight) )
    gs = gridspec.GridSpec( rows, 2, wspace=wspace, hspace=hspace, width_ratios=(ratio_left,1-ratio_left) )

    # generate axis array:
    ax = np.empty((rows,2), dtype=object)
    ax[rows-1,0] = plt.subplot( gs[2*rows-2] ) # first the bottom plot to share the axis 
    shareyax = None if not sharey else ax[rows-1,0]
    ax[rows-1,1] = plt.subplot( gs[2*rows-1], sharey=shareyax ) # first the bottom plot to share the axis
    for axis in range(0,rows-1): # then the other axes
        ax[axis,0] = plt.subplot( gs[axis*2]  , sharex = ax[rows-1,0] ) # share axis
        shareyax = None if not sharey else ax[axis,0]
        ax[axis,1] = plt.subplot( gs[axis*2+1], sharex = ax[rows-1,1], sharey=shareyax ) # share axis
        ax[axis,0].tick_params( 'x', labelbottom = False ) # hide ticks 
        ax[axis,1].tick_params( 'x', labelbottom = False ) # hide ticks

    # set ticks:
    if ticks == "everywhere":
        for axis in range(0,rows):
            ax[axis,1].yaxis.tick_right() # ticks to the right
            ax[axis,0].tick_params( right = True )
            ax[axis,1].tick_params( left = True )
            ax[axis,0].tick_params( right = True, which = 'minor' )
            ax[axis,1].tick_params( left = True, which = 'minor'  )
    elif ticks == "half":
        for axis in range(0,rows):
            ax[axis,1].yaxis.tick_right() # ticks to the right
    else:
        raise Exception("tick setting %s unknown"%ticks)

    # set labels:
    for axis in range(0,rows):
        ax[axis,1].yaxis.set_label_position("right") # labels to the right

    return fig, ax
