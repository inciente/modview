# Functions to find index mapping between data objects of different types.
import xarray as xr; import numpy as np; 

def var_at_iso(varq, varl, values):
    # Take in two xr.DataArrays varq, varl, and map the values in varq onto 
    # lines/surfaces where varl == values[kk]. 
    if varq.dims == varl.dims:
        print('var_at_iso has not been written yet')
        pass
    else:
        print('var_at_iso has not been written yet')    

def find_in_arr(arr, value):
    # Return the index of item in arr that is closest to value.
    delta = abs(arr-value);
    return np.argmin(delta) 
    
def find_axis(var, vec):
    ''' Return axis along which var has the same size as vec
    '''
    # dimensions of var with same length as vec
    check = [var.shape[kk] == len(vec) for kk in range(len(var.shape))]; 
    axis = [i for i, x in enumerate(check) if x]; # store dimension numbers 
    if len(axis)>1:
        print('dimensions agree along more than one axis')
        return
    else:
        return axis[0]
    
def put_in_bins(data, coords, bins):
    ''' data is a 1D array that represents nd data, 
    coords is a list with n arrays, each with coordinates entries in data
    bins is also a list with n arrays
    fn
    '''
    if np.size(coords) != np.size(bins):
        print('Coords and bins must be lists with same number of elements')
        break 
    in_bins = []; 
    for sdim in range(len(coords)):
        check = np.digitize(coords[sdim], bins[sdim]); # which entry in which bin
        
    # Take the average of values in obj_from between     
    axis = find_axis(obj_from, coord_from); 
    print('avg_within_bins has not been written yet')

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

def match_latlon(obj, to_obj):
    newobj = obj.sel(lon=to_obj.lon);
    newobj = obj.sel(lat=to_obj.lat); 
    return newobj

def middlepoint(arr):
    # Return the item in the middle of a n-dimensional np array.
    shp = np.floor(np.array(arr.shape)/2); 
    shp = [[int(kk)] for kk in shp]; # list of lists for fancy indexing
    val = arr[shp];
    return val
