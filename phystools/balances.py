import numpy as np; import xarray as xr; import pandas as pd; 
import sympy as sym; import gsw, datetime
from abc import ABC, abstractmethod
from sympy.utilities.lambdify import lambdify, implemented_function
from . import mapper

g = 9.81; rho_0 = 1024; C_p = 3.9e3; 

def inertial_freq(lat):
    f = 4*np.pi*np.sin(lat/180*np.pi)/24/3600
    return f

def beta_f(lat):
    beta = 4*np.pi*np.cos(lat/180*np.pi)/24/3600;
    beta = beta/6.371e6; # divide by Earth's radius
    return beta

def hydrostasy(rho, zlevels, ssh=0, z_req=np.nan, p_axis=0, p_ref=0):
    # Taake in n-dimensional object rho and integrate pressure downwards
    # ssh represents local values at locations of density data
    # currently only works for 1d objects of density
    dz = np.diff(zlevels); 
    pressure = np.empty(rho.shape); 
    pressure[:] = np.nan; 
    pressure[0] = (ssh + zlevels[0])*g*min(rho); # pressure at shallowest level
    pressure[1:] = [dz[i]*g*(rho[i-1]+rho[i])/2 for i in range(1,len(pressure))]; # p diff btween levels
    pressure = np.cumsum(pressure); 
    return pressure

def wkb_stretch(z,N2,MLD):
    belowML = z>MLD;
    N = np.sqrt(N2);
    sratio = np.nanmean(N[belowML])/N;
    sratio = pd.DataFrame(sratio); # to get rolling mean
    sratio = sratio.rolling(10, min_periods=1).mean().to_numpy();

    vel_factor = np.sqrt(sratio); # velocity scaling

    tot_depth = z[~belowML]; # stretch vertical coordinate
    z_stretch = np.expand_dims( np.append(0, np.diff(z)), axis=1);
    z_stretch = np.cumsum( z_stretch[belowML]/sratio[belowML])+MLD; # vertical integration
    z_stretch = np.append(tot_depth, z_stretch);
    return vel_factor, stretched

def lat_lon_deriv(variable, location):
    ''' Put in an xr_object and a coordinate dictionary (location) and calculate
    derivatives on x and y for variable. Check past implementations to know exactly how
    location must be defined and whether or not it can be obtained directly from variable'''
    var_x = []; #xr.DataArray(data=np.nan, coords=location); 
    var_y = []; #xr.DataArray(data=np.nan, coords=location); 
    if 'lat' in variable.dims:
        var_y = variable.differentiate(coord='lat')/110e3;
    if 'lon' in variable.dims:
        var_x = variable.differentiate(coord='lon')/110e3 \
                / np.cos(np.radians(location['lat'])); 
    return var_x, var_y

class geostrophic_field:
    def __init__(self,CT,SAL,P,ssh='none',dims='both'):
        # all inputs must be xr.dataarrays. ssh may have different coordinates
        self.CT = CT
        self.SAL = SAL; 
        self.P = P; 
        self.ssh = ssh; 
        # now calculate some basic quantities
        self.p_axis = self.find_axis(self.SAL, self.P); # might want to move this function to mapper
        self.density = gsw.density.rho(self.SAL, self.CT, self.P); 
        self.psi = self.streamfunction(); 

    def streamfunction(self):
        p_mat = self.upsize_p(self.CT); 
        # calculate dynamic height
        DH = gsw.geo_strf_dyn_height(self.SAL.values, self.CT.values, p_mat, axis=self.p_axis); 
        DH = xr.DataArray( data=DH, dims=self.CT.dims, coords=self.CT.coords); 
        return DH

    def geost_flow(self,location,edge=1,use='psi',point=False):
        # calculate u for the entire field using ssh and density (need to copy from phystools.py)
        pass

    # MOVE THIS ONE TO MAPPER
    def var_around(self, xr_obj, location, edge=1):
        subsel = mapper.block_around(xr_obj, location, edge=edge); 
        return subsel

