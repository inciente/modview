# Functions to find index mapping between data objects of different types.
import xarray as xr; import numpy as np; 

def var_at_iso(varq, varl, values):
    # Take in two xr.DataArrays varq, varl, and map the values in varq onto 
    # lines/surfaces where varl == values[kk]. 
    if varq.dims == varl.dims:
        pass
    else:
        print('fake')    

def find_in_arr(arr, value, method='min'):
    # Return the index of item in arr that is closest to value.
    delta = abs(arr-value);
    return np.argmin(delta) 
   
def upsize_repeat(block, vec):
    # block and vec must have one common dimension
    # return a block that repeats vec such that block, return have same shape
    match_axis = find_axis(block,vec); 
    new_block = np.ones(block.shape); 
    # Dimensions missing in vec
    ax_to_create = [kk for kk in range(len(block.shape)) if kk!=match_axis]; 
    p_col = np.expand_dims( vec, tuple(ax_to_create)); # create new dimensions
    new_block = new_block*p_col; 
    return new_block

def find_axis(var, vec):
    # Return axis along which var has the same size as vec
    check = [var.shape[kk] == len(vec) for kk in range(len(var.shape))]; 
    axis = [i for i, x in enumerate(check) if x]; 
    if len(axis)>1:
        print('dimensions agree along more than one axis')
        return
    else:
        return axis[0]
    
def avg_within_bins(obj_from, coord_from, coord_to):
    # Take the average of values in obj_from between     
    axis = find_axis(obj_from, coord_from); 
    

def block_around(xr_obj, location, edge=1, time=False):
    # Take in an xr.dataarray (xr_obj), and extract a block of data around 
    # the location stated in dictionary coords. 
    if tuple([key for key in location]) != xr_obj.dims:
        print('dimensions dont match, bro')
        return
    # Add if condition to select subset of xr_obj for corresponding time.    
    indices = []; dim_num = 0; 
    new_coords = dict();  
    for key in location:
        if isinstance(edge,int): # how big is the block in each direction?
            d_edge = np.arange(-edge, edge+0.1,1); 
        elif isinstance(edge,list):
            d_edge = np.arange(-edge[dim_num], edge[dim_num]+0.1, 1); 
        dim_num += 1
        coord_ax = xr_obj[key].values;
        central = find_in_arr(coord_ax, location[key]); # point closest to requested location
        get_items = [int(central+dl) for dl in d_edge]; # add edge indices
        indices.append( get_items );
        new_coords[key] = coord_ax[get_items];    
    ind2get = np.meshgrid(*indices, indexing='ij'); # block-like set of indices
    subsel = xr_obj.values[ind2get];
    subsel = xr.DataArray( data=subsel, dims=xr_obj.dims, coords=new_coords); 
    return subsel


def standardize_coords(obj):
    ''' Now that standard coordinate names have been set in coord_axes(), 
    create a function that will rename dimensions in any given xr using the 
    axes dictionary that comes out of coord_axes. That way the dictionary 
    does not have to be carried over within external functions and classs. '''
    axes = coord_axes(obj); 
    renamer = dict(); 
    for st_coord in axes.keys():
        renamer[axes[st_coord][0]] = st_coord;
    obj.rename(renamer);
    return obj

def coord_axes(obj):
    ''' Read the names and values of axes/dims in xr obj and
    return a dictionary that identifies dimensions with name and index.
    This will help treat all geographical objects using the standardized 
    coordinates: longitude, latitude, pressure, time (x,y,z,t).'''
    ND = len(obj.shape); # number of dimensions
    data_dims = obj.coords.dims; # tuple of strings
    axes = dict(); 
    # -----------------------------------------------------
    for dim in data_dims:
        # save coord as string, values (xr), and dimension axis in obj
        translator = (dim, obj[dim], find_axis(obj, obj[dim])); 
        # save in corresponding dictionary key
        if dim in ['lon','longitude','LONGITUDE','LON']:
            axes['longitude'] = translator; 
        elif dim in ['lat','latitude','LATITUDE','LAT']:
            axes['latitude'] = translator; 
        elif dim in ['p','pressure','PRESSURE','depth']:
            axes['pressure'] = translator; 
        elif dim in ['time','date','t']:
            axes['time'] = translator
        else: 
            print('Dimension ' + dim + ' was not found in object');
    return axes

def match_latlon(obj, to_obj):
    newobj = obj.sel(lon=to_obj.lon, method='linear');
    newobj = obj.sel(lat=to_obj.lat, method='linear'); 
    return newobj

def middlepoint(arr):
    # Return the item in the middle of a n-dimensional np array.
    shp = np.floor(np.array(arr.shape)/2); 
    shp = [[int(kk)] for kk in shp]; # list of lists for fancy indexing
    val = arr[shp];
    return val

def coords_to_xyzt(xr_obj):
    ''' Take in an xr object and return a dictionary of coordinates
    that can replace the original object's coords. '''
    use_dims = xr_obj.coords.dims; 
    new_coords = dict(); 
    new_coords[x] = xr_obj['longitude'] * 110e3 \
                 / np.cos( np.nanmean(xr_obj['latitude'])/180*np.pi);
    new_coords[y] = xr_obj['latitude'] * 110e3; 
    if 'pressure' in use_dims:
        new_coords[z] = -xr_obj['pressure'];
    if 'time' in use_dims:
        # call timetools
        pass
    return new_coords















