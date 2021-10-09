import numpy as np; import xarray as xr; import pandas as pd; 
import sympy as sym;
from sympy.utilities.lambdify import lambdify, implemented_function

def lon_to_x(lon, lat):
    x = lon*110e3*np.cos(np.mean(lat)); x = x - np.mean(x);
    return x

def lat_to_y(lat): 
    y = lat*110e3; y = y - np.mean(y);
    return y
    
def inertial_freq(lat):
    f = 4*np.pi*np.sin(lat/180*np.pi)/24/3600;
    return f
   
def geostrophy(lat, lon, depth, density, ssh='none'):
    # Calculate geostrophic shear/currents from density along a transect.
    # either lat or lon  need to be sequence with a single element
    if len(lat)>1 and len(lon)>1:
        print('Coordinates dont correspond to a transect')
    elif len(lat) == 1 and len(lon) > 1: 
        solve_for = 'v';
    elif len(lon) == 1 and len(lat) > 1:
        solve_for = 'u';
    f = 4*np.pi*np.sin(np.mean(lat))/24/3600;
    x = lon_to_x(lon,lat); 
    y = lat_to_y(lat);    
    
    if density.shape[0] == len(depth):
        drho = np.gradient( density, axis=1); 
    elif density.shape[1] == len(depth):
        drho = np.gradient( density, axis=0); 
    else: 
        print('Depth and density dimensions dont match')
    
    if solve_for == 'u':
        drho = drho / np.gradient(y); 
    elif solve_for == 'v':
        drho = drho / np.gradient(x); 

class internal_wave:
    def __init__(self, lat, N2=1e-5 ):
        self._fk = dict();  # vector [freq, kx, ky, kz].
        self._Nsquared = N2;
        self.lat = lat;
    
    @property
    def Nsquared(self):
        print('getting N2');
        return self._Nsquared
    
    @Nsquared.setter
    def Nsquared(self, n2):
        # Set value to argument n2
        self._Nsquared = n2;
    
    @property
    def fk(self):
        return self._fk
    
    def fk_sym(self):
        wave_properties = [(sym.symbols(key),value) for key, value in self._fk.items()]
        return wave_properties
    
    @fk.setter
    def fk(self, wave_properties):
        # Run dis_rel using fk_vec input. 
        test_func = self.dispersion_relation()
        
        test_input = []; solve_for = []; var_solved = [];
        for key, value in wave_properties.items():
            if np.isnan(value):
                    var_solved.append(key);
                    solve_for.append(sym.symbols(key)); # values missing
            else: 
                test_input.append( (sym.symbols(key), value) ) # values known
        try_func = test_func.subs(test_input);
        solution = sym.solve(try_func, solve_for)
        
        if solution[0].is_real:
            wave_properties[var_solved[0]] = solution[0];
            self._fk = wave_properties; 
        else:
            print('Incompatible values: imaginary wavenumber');
            pass

    def dispersion_relation(self, numeric=False):
        if numeric:
            
            
        # Symbolic representation of dispersion relation
        omega, kx, ky, kz = sym.symbols('omega kx ky kz')
        norm = kx**2 + ky**2 + kz**2;
        relation = - omega**2 + self.Nsquared*(kx**2 + ky**2)/norm \
                    + inertial_freq(self.lat)**2 * kz**2 / norm; 
        return relation
    
    def group_vel(self, direction='z'):   
        kz = sym.symbols('kz');
        cg_z = sym.diff(self.dispersion_relation(),kz)
        cg_val = sym.lambdify(self.fk_sym(),cg_z)
        return cg_val

        
