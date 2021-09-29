import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import scipy.interpolate as interpolate
import pandas as pd

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
    
def grid_xyz(df,y,y_query,dt_query,limits='none', log=False):
	# Take in a pandas dataframe with values (z) and index (x). 
	# y is 1d vector or pd.df with same shape as df. 
	# dt_query is a string (argument of resample)
    if isinstance(limits,dict):
        df = df[limits['t0']:limits['t1']]; 
        y = y[limits['t0']:limits['t1']]; 
    df = df.resample(dt_query).mean();
    if log:
        df = np.log10(df);  
    y = y.resample(dt_query).mean(); 
    gridded = np.empty([df.shape[0],len(y_query)]);
    gridded[:] = np.nan;
    for tt in range(df.shape[0]):
        locations = y.values[tt,:]; 
        valz = df.values[tt,:]; 
        finterp = interpolate.interp1d( locations, valz)
        gridded[tt,:] = finterp(y_query);
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
