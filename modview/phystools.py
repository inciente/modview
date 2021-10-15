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
    
    def fk_vals(self): 
        wave_properties = self.fk;
        kx = wave_properties['kx']; ky = wave_properties['ky']; 
        kz = wave_properties['kz']; omega = wave_properties['omega'];
        return omega, kx, ky, kz
    
    @fk.setter
    def fk(self, wave_properties):
        # Run dispersion_relation using wave_properties input. 
        test_func = self.dispersion_relation()
        
        test_input = []; solve_for = []; var_solved = [];
        for key, value in wave_properties.items():
            if np.isnan(value):
                    var_solved.append(key);
                    solve_for.append(sym.symbols(key)); # values missing
            else: 
                test_input.append( (sym.symbols(key), value) ) # values known
        try_func = test_func.subs(test_input);        
        if solve_for == []: # if all values are known
		    # evaluate try_func
            if try_func == 0:
                self._fk = wave_properties;
            else:
                print('wave properties not compliant with dispersion relationship')
                pass
        else:	
            solution = sym.solve(try_func, solve_for) # use known values to infer missing values
        
            if solution[0].is_real:
                wave_properties[var_solved[0]] = solution[0];
                length = 2*np.pi/solution[0]; 
                print('Missing parameter ' + var_solved[0] + '= 2*pi/' \
                + '{:5.2f}'.format(length) + ' added successfully')
                self._fk = wave_properties; # update wave properties
            else:
                print('Incompatible values: imaginary wavenumber');
                pass

    def dispersion_relation(self, numeric=False):
        if numeric: # evaluate terms in the dispersion relation
            wave_properties = self.fk;
            kh_squared = wave_properties['kx']**2 + wave_properties['ky']**2;
            norm_squared = kh_squared + wave_properties['kz']**2; 
            om2 = wave_properties['omega']**2; 
            t1 = self.Nsquared * kh_squared / norm_squared;
            t2 = inertial_freq(self.lat)**2 * wave_properties['kz']**2 / norm_squared;
            
            relation = [[om2, t1, t2],['omega^2','N^2 kh^2/k^2','f^2 kz^2/k^2']]; 
            
        else:    # Symbolic representation of dispersion relation
            omega, kx, ky, kz = sym.symbols('omega kx ky kz')
            kh_squared = kx**2 + ky**2;
            norm_squared = kh_squared + kz**2;
            relation = - omega**2 + self.Nsquared*kh_squared/norm_squared \
                    + inertial_freq(self.lat)**2 * kz**2 / norm_squared; 
        return relation
    
    def group_vel(self, direction=['z'], numeric=False, kunze=False):   
        cg_vec = [np.nan, np.nan, np.nan];
        kx,ky,kz = sym.symbols('kx ky kz'); # symbols for sympy functions
        disp_rel = self.dispersion_relation();
        
        if 'x' in direction: # compute only necessary derivatives
            cg_vec[0] = sym.diff(disp_rel,kx)
        if 'y' in direction:
            cg_vec[1] = sym.diff(disp_rel,ky)
        if 'z' in direction:
            cg_vec[2] = sym.diff(disp_rel,kz)
        
        if numeric: # substitute wave properties into 
            wave_input = self.fk_sym();
            if 'x' in direction:
                cg_vec[0] = cg_vec[0].subs(wave_input);
            if 'y' in direction: 
                cg_vec[1] = cg_vec[1].subs(wave_input);
            if 'z' in direction: 
                cg_vec[2] = cg_vec[2].subs(wave_input); 
        return cg_vec

    def is_niw(self, threshold=0.2):
        inertial = inertial_freq(self.lat); 
        freq = self.fk['omega']; 
        ratio = freq/inertial; 
        if ratio > (1-threshold) and ratio < (1+threshold):
            return True
        else: 
            return False
        
    def setup_kunze(self, dudy, dvdx, dudz=0, dvdz=0, dudx=0, dvdy=0):    
        # input necessary: dU/dz, dV/dz, dU/dy, dV/dx
        if self.is_niw():
            f_loc = inertial_freq(self.lat); 
            vort = dvdx - dudy; 
            f_eff = f_loc + 0.5*vort;
            omega, kx, ky, kz = self.fk_vals(); 
            kh_squared = kx**2 + ky**2; 
            omega_0 = f_eff + ( # term (i)
            # term (ii)
                 + self.Nsquared*kh_squared/(2*f_loc*kz**2)) + (
            # term (iii)
                  1/kz*( dudz*ky - dvdz*kx)) + ( 
            # term (iv)      
                  1j*( vort/2/f_loc/kz*(dudz*kx + dvdz*ky)  
            # term (v)
                    + self.Nsquared/(2*f_loc**2*kz**2)*( 
                         dudx*ky**2 - (dudy + dvdx)*kx*ky + 
                         dvdy*kx**2 ) ))  
            print(omega_0)
        else:
            print('wave is not near-inertial')
