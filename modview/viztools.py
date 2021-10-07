import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import scipy.interpolate as interpolate
import pandas as pd

def dateticks(axlist, axind, dts):
    # Major ticks.
    fmt_big = mpl.dates.DayLocator(interval=dts[1])
    axlist[axind].xaxis.set_major_locator(fmt_big)
    # Minor ticks.
    fmt_small = mpl.dates.DayLocator(interval=dts[0])
    axlist[axind].xaxis.set_minor_locator(fmt_small)
    # FOrmat of text
    axlist[axind].xaxis.set_major_formatter(mpl.dates.DateFormatter('%b %-d'))

def dateticks2(ax_obj, dts):
    # Major ticks.
    fmt_big = mpl.dates.DayLocator(interval=dts[1])
    ax_obj.xaxis.set_major_locator(fmt_big)
    # Minor ticks.
    fmt_small = mpl.dates.DayLocator(interval=dts[0])
    ax_obj.xaxis.set_minor_locator(fmt_small)
    # FOrmat of text
    ax_obj.xaxis.set_major_formatter(mpl.dates.DateFormatter('%b %-d'))

def scatter_xyz(x,y,color_var,y_offset='none',m_axis=1,vlims='auto'):
	# Input a matrix (x). It is assumed that each column (m_axis=1) is a 
	# different set of measurements to be plotted.
	# y is a vector or matrix such that (x[:,0], y[0]+y_offset) or 
	# (x[:,0], y[:,0]) are to be plotted.
    if x.shape == color_var.shape:
        x2plot = np.reshape(x,-1);
        c2plot = np.reshape( color_var,-1); 
    else:
        print('x and c not same shape');
    if x.shape == y.shape:
        y2plot = np.reshape(y,-1); 
    else:
        # This is where we add y_offset. y should be (cols(x),1) and 
        # y_offset (1, rows(x))
        newy = y + y_offset; 
        y2plot = np.reshape(newy,-1); 
    sf = plt.scatter(x2plot, y2plot, c=c2plot, vmin=vlims[0], vmax=vlims[1]); 
    return sf
    
    
def grid_xyz(df,y,y_query,dt_query,limits='none', log=False, data_format='xr'):
	# Take in a pandas dataframe with values (z) and index (x). 
	# y is 1d vector or pd.df with same shape as df. 
	# dt_query is a string (argument of resample)
    if data_format=='pd':
        if isinstance(limits,dict):
            df = df[limits['t0']:limits['t1']]; 
            y = y[limits['t0']:limits['t1']]; 
        df = df.resample(dt_query).mean();
        y = y.resample(dt_query).mean(); #movement of chipods
        if log:
            df = np.log10(df);  
        
        gridded = np.empty([df.shape[0],len(y_query)]);
        gridded[:] = np.nan; #interpolate onto this empty matrix
        for tt in range(df.shape[0]): # for each timestep
            locations = y.values[tt,:]; # read chipod depths
            valz = df.values[tt,:]; # measurement
        
            finterp = interpolate.interp1d( locations, valz, bounds_error='none');
            # Check for query inside observation range
            gridded[tt,:] = finterp(y_query);
        if tt == 40:
            plt.scatter(valz, locations,s=30)
            plt.plot(finterp(y_query), y_query)  
            print(valz)
            print(locations)
    gridded = pd.DataFrame(data=gridded, index=df.index, columns=y_query);      
    return gridded
    
class MplColorHelper:

  def __init__(self, cmap_name, start_val, stop_val):
    self.cmap_name = cmap_name
    self.cmap = plt.get_cmap(cmap_name)
    self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
    self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

  def get_rgb(self, val):
    return self.scalarMap.to_rgba(val)
