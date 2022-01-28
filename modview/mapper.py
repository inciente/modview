# Functions to find index mapping between data objects of different types.
import xarray as xr; import numpy as np; 
from gsw import distance 
import datetime
from abc import ABC, abstractmethod

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


def standardize_coords(obj,array=None):
    ''' Now that standard coordinate names have been set in coord_axes(), 
    create a function that will rename dimensions in any given xr using the 
    axes dictionary that comes out of coord_axes. That way the dictionary 
    does not have to be carried over within external functions and classs. '''
    if isinstance(obj, xr.DataArray):
        axes = coord_axes(obj);
    elif isinstance(obj, xr.Dataset):
        axes = coord_axes( obj[array] );
    renamer = dict(); 
    for st_coord in axes.keys():
        renamer[axes[st_coord][0]] = st_coord;
    return obj.rename(renamer)


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
        if dim in ['lon','longitude','LONGITUDE','LON','X']:
            axes['longitude'] = translator; 
        elif dim in ['lat','latitude','LATITUDE','LAT','Y']:
            axes['latitude'] = translator; 
        elif dim in ['p','pressure','PRESSURE','depth','P']:
            axes['pressure'] = translator; 
        elif dim in ['time','date','t','TIME','T']:
            axes['time'] = translator
        else: 
            print('Dimension ' + dim + ' was not found in object');
    return axes

def mean_line(array,axis):
    ''' Take in an array and extract a 1D slice of all elements along
    a given axis. This was written to turn 3D and 4D pressure data that
    comes out of gsw.Nsquared() into usable 1D arrays.'''
    ndims = len(array.shape); 
    for jj in np.flip(range(ndims)):
        if jj != axis:
            array = np.nanmean(array, axis=jj);
    array = np.squeeze(array); 
    return array

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


class location:
    def __init__(self, longitude=np.nan, latitude=np.nan, \
            pressure=np.nan, time=np.nan):
        # Create brand new position.
        self.loc = np.array([longitude, latitude, pressure, time]); 
        #self.longitude = longitude; 
        self.latitude = self.loc[1]; 
        #self.pressure = pressure; 
        #self.time = time;     
        
    def distance_to(self, other_location):
        obj_list = [self, other_location];
        dl = distance([kk.loc[0] for kk in obj_list],
                [kk.loc[1] for kk in obj_list], 
                [kk.loc[2] for kk in obj_list]); # in meters 
        return dl

    def displace(self, **kwargs):
        ''' Returns a new location by adding displacements to each 
        coordinate in self.loc. Displacements can be given in meters (dx,dy)
        or degrees (dlon, dlat).'''
        disp = np.array([0 if ~np.isnan(kk) else np.nan for kk in self.loc]); 
        midlat = self.loc[1];
        if 'dy' in kwargs:
            disp[1] = kwargs['dy']/110e3;
        elif 'dlat' in kwargs:
            disp[1] = kwargs['dlat'];
            midlat = self.loc[1] + dlat/2
        if 'dx' in kwargs:
            disp[0] = kwargs['dx']/110e3/np.cos(midlat/180*np.pi);
        elif 'dlon' in kwargs:
            disp[0] = kwargs['dlon']; 
        if 'dp' in kwargs:
            disp[2] = kwargs['dp'];
        if 'dt' in kwargs:
            disp[3] = kwargs['dt'];
        new_loc = self.loc + disp;
        new_loc = location( longitude=new_loc[0], latitude=new_loc[1], 
                pressure=new_loc[2], time=new_loc[3]); 
        return new_loc

    def eval_xr(self,xr_obj, method='linear', dlon=0, dlat=0, dp=0):
        if self.loc[2] is not None and \
                'pressure' in xr_obj.coords.dims:
            xr_obj = xr_obj.interp(pressure=self.loc[2]+dp, method=method);
        xr_obj = xr_obj.interp(longitude=self.loc[0]+dlon, method=method);
        xr_obj = xr_obj.interp(latitude=self.loc[1]+dlat, method=method); 
        return xr_obj.values

    def gradient_xr(self, xr_obj, dlat=None, dlon=None, dp=None):
        ''' Calculate gradient of variable in an xr.Dataarray
        at this location. Depending on which arg is not none, we
        will calculate gradient in a different direction.'''
        grad = [None, None, None]; 
        if dlon is not None:
            f1 = self.eval_xr(xr_obj, dlon=-dlon/2); 
            f2 = self.eval_xr(xr_obj, dlon=dlon/2); 
            grad[0] = (f2-f1)/(dlon*110e3*np.cos(self.latitude))
        if dlat is not None:
            f1 = self.eval_xr(xr_obj, dlat=-dlat/2);
            f2 = self.eval_xr(xr_obj, dlat=dlat/2); 
            grad[1] = (f2-f1)/(dlat*110e3)
        if dp is not None:
            f1 = self.eval_xr(xr_obj, dp=-dp/2); 
            f2 = self.eval_xr(xr_obj, dp=dp/2); 
            grad[2] = -(f2-f1)/dp # relative to z (increase upward)
        return grad

''' DEFINE FACTORY FOR LOCATIONS AND CALL FACTORIES FROM WITHIN TRACK 
that way, track can create new locations without having to know anything about
how those locations are being generated.

Abstract base classes needed: 

'''

class sampling_scheme(ABC):
    ''' Subclasses may include: static (constant set of locations), 
    track(s), where locations changes over time
    '''
    
    def __init__(self, **kwargs):
        self.elements = dict(); 

    def get_element(self,el_id):
        print('getting element from sample')

    def create_location(self, **kwargs) -> location:
        nu_loc = location(self, kwargs); 
        return nu_loc

   # def assign_loc_to_element(self,element,location):





class track:
    def __init__(self, list_of_locations=[]):
        self.path = list_of_locations; 
        self.dims = []; 

    def add_location(self,loc):
        if isinstance(loc,location):
            self.path.append(location); 
        else: 
            print('Warning: items must be instances of the location class');

    def to_xyzt(self):
        ''' Return arrays with the values of longitude, latitude, etc. 
        To simplify plotting time series '''
        pos_dict = {'longitude':[], 'latitude':[],
                'pressure':[], 'time':[]};
        index = 0
        for loc in self.path:
            pos_dict['latitude'].append(loc.loc[1]);
            pos_dict['longitude'].append(loc.loc[0]); 
            pos_dict['pressure'].append(loc.loc[2]); 
            pos_dict['time'].append(loc.loc[3]);
        for kk in pos_dict.keys():
            pos_dict[kk] = np.array(pos_dict[kk]); 
        return pos_dict











