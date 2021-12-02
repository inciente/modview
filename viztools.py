import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import scipy.interpolate as interpolate
import pandas as pd
import datetime

""" Use this class to create multi-panel visualizations.
Inputs are two dictionaries. 
graphic_dict must include: ['figsize','widths','heights','panels'] to specify
the general format of the plot.
""" 

class panel_plot:
    def __init__(self,graphic_dict, data_dict='none'):
        self.data = data_dict;
        self.graphic = graphic_dict;
        self.axes = []; 
        self.fig = plt.figure(figsize=graphic_dict['figsize']); 
        self.make_gs();
        self.make_axes();
        self.write_labels(); 
               
    # What do I want from this class?
    # Create grid using specifications in graphic_dict
    # Find a way to relate data to specific panels so it can all be added where it needs to go.
    # Create grid spec using figure information
    
    def make_gs(self): 
        nx = len(self.graphic['widths']); 
        ny = len(self.graphic['heights']); 
        gs = self.fig.add_gridspec(ncols=nx, nrows=ny, 
                    width_ratios=self.graphic['widths'], 
                    height_ratios=self.graphic['heights']); 
        self.graphic['gs'] = gs; 

    def make_axes(self):
        panlist = self.graphic['panels']; 
        for kk in panlist:
            ax_here = self.fig.add_subplot( self.graphic['gs'][kk[0],kk[1]]);
            self.axes.append(ax_here)
    
    def write_labels(self):
        jj=0
        for ax in self.axes:
            ax.set_ylabel(self.graphic['labels'][jj])
            jj+=1
        
        
    def paint_panel(self,axind,datloc, **kwargs):
        # Use datloc to locate data in data_dict to paint it on graphic_dict['axes'][axind]. 
        # Use **kwargs to specify the format of the graphic  
        pass
    
    def date_axis(self, axind, limits, ticks, which='x'):
        # Set equal axis limits for all panels in list axind. 
        for jj in range(len(axind)):
            self.axes[jj].set_xlim([pd.to_datetime(limits['t0']), pd.to_datetime(limits['t1'])])
            dateticks2(self.axes[jj], ticks)
       

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
    
def cut_var(obs, dsource, variable):
    if dsource =='ADCP':
        var_dat = obs.vars[variable];
        var_dat = var_dat.sel(z=slice(limits['z0'],limits['z1']));
        z_here = var_dat['z'].values;
    elif dsource == 'chipod':
        var_dat = obs.vars['chipod'][variable]
        #var_dat = var_dat.sel(z=slice(limits['z0'],limits['z1']));
        z_here = obs.vars['chipod']['z'];
    elif dsource == 'Temperature':
        var_dat = obs.vars[dsource][variable]; 
        z_here = var_dat['depth']; 
    var_dat = var_dat.sel(time=slice(limits['t0'],limits['t1'])); # apply limits
    return var_dat
        
    
def plot_period(obs, dsource,variable, axes, limits,resamp='none', log=False, viztype='pcolor',**kwargs):
    # obs is an instance of loader.assemble
    # variables is dict() to access variable to plot
    # axes is where the plot will be made. 
    # limits is the typical dict with t0,t1,z0,z1
    var_dat = cut_var(obs, dsource, variable, limits); 
    
    if resamp == 'none':
        pass
    else:
        var_dat = var_dat.resample(time=resamp).mean(); 
    if log: # get data in plotting form
        var_vals = np.log10(var_dat.values);
    else: 
        var_vals = var_dat.values;
    if dsource=='chipod':
        var_vals = var_vals.transpose();
        
    # Prepare visualization format 
    viz_args = {'cmap':'viridis','vmin':-1,'vmax':1,'shading':'flat','levels':'none'};
    for arg in viz_args.keys():
        if arg in kwargs:
            viz_args[arg] = kwargs.get(arg);
            
    tvec = pd.to_datetime(var_dat.time.values)
    if viztype=='contourf':
        del viz_args['vmin'], viz_args['vmax']
        panel=axes.contourf(tvec, z_here,var_vals, **viz_args)
    elif viztype=='pcolor':
        del viz_args['levels']
        panel = axes.pcolormesh(tvec, z_here, var_vals, **viz_args)
    elif viztype=='contour':
        del viz_args['vmin'], viz_args['vmax']
        panel=axes.contour(tvec,z_here, var_vals, **viz_args)
    return panel

    
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
