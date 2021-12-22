''' Create meta classes for waves and ray tracing methods that couple waves to classes of
medium fields (beta plane, geostrophic flow, etc.). Methods will rely on both symbolic represntations
of wave properties and numerical representations of media. Use modview.mapper '''


import pandas as pd; import numpy as np; import datetime, math; 
import sympy as sym; 
from abc import ABC, abstractmethod
import sys, importlib
sys.path.append('/mnt/sda1/PhysOc/modview/')
import mapper

# Begin by declaring some useful functions
def inertial_freq(lat):
    f = 4*np.pi*np.sin(lat/180*np.pi)/24/3600
    return f

class wave(ABC):
    ''' Meta class to represent different waves in geophysical fluid dynamics.
    Main attributes are:
        1. Wave vector - dictionary including frequency and wavenumbers
        2. Medium - dictionary including location and environmental variables in dispersion relationship
        3. Dispersion relationship - defind symbollically within subclasses.
    '''
    def __init__(self):
        pass

    @property 
    def wave_vector(self): 
        return self._wave_vector
    
    @wave_vector.setter
    def wave_vector(self,value):
        test, fk = self.validate_wvnum(fk=value) # does input comply with disp rel?
        if test:
            self._wave_vector = fk; # set property
            print(fk)
        else: 
            print('Warning: wave vector does not satisfy dispersion relationship!')
            sys.exit()
           
    @property
    @abstractmethod
    def medium(self):
        return self._medium
    
    @medium.setter
    def medium(self,value):
        # Extract all variables (f, beta, N^2, depth, etc.) that characterize 
        # the propagation medium and influence wave properties/behavior
        pass

    # learn more about abstract properties.
    @abstractmethod
    def dispersion_rel(self,numeric,**kwargs):
        pass 
    
    @abstractmethod
    def validate_wvnum(self, fk):
        ''' Test whether a proposed wave vector dictionary fk 
        complies with the dispersion relationship. 
        Abstract for now, but might be implemented in metaclass in the future.'''
        pass
        
    def sym_vector(self):
        ''' Return a list of symbols representing the frequency and 
        wavenumbers of a given wave.'''
        fk_sym = [(sym.symbols(key),value) for key, 
                       value in self.wave_vector.items()]
        return fk_sym
    
    def check_type(self,**kwargs):
        if 'fk' in kwargs.keys():
            fk = kwargs['fk'];
        else:
            fk = self.wave_vector;
        disp_terms = self.dispersion_rel(numeric=True,fk=fk) # evaluate dispersion relation
        test = math.isclose(disp_terms[0][0], disp_terms[0][1], 
                               rel_tol=1e-6) # compare left and right sides
        return test

class internal(wave):
    def __init__(self):
        self._wave_vector = dict()
        self._medium = dict()
    # Create symbolic representation of ...
    omega, kx, ky, kz = sym.symbols('omega kx ky kz'); # wave properties
    N, f = sym.symbols('N f'); # environmental variables
    norm2 = kx**2 + ky**2 + kz**2;  
    lhs = omega; 
    rhs = sym.sqrt( N**2 * (kx**2 + ky**2)/norm2 + f**2 * kz**2 / norm2 ); 

    c_g = [sym.diff(rhs, kx), sym.diff(rhs, ky), sym.diff(rhs, kz)];

    @property 
    def medium(self):
        return self._medium
    
    @medium.setter
    def medium(self,val_dict):
        self._medium = val_dict; # include lat, f, others

    def validate_wvnum(self,fk):
        disp_rel = -self.lhs + self.rhs; 
        solve_for = []; var_solved = [];
        test_input = [(self.N, np.sqrt(self.medium['N2'])), 
                      (self.f, inertial_freq(self.medium['lat'])) ];
        # separate known from unknown values
        for key, value in fk.items():
            if np.isnan(value):
                var_solved.append(key);
                solve_for.append(sym.symbols(key)); # missing value
            else:
                test_input.append( (sym.symbols(key), value) ); # known values
        
        try_func = disp_rel.subs(test_input); # substitute known values
        if solve_for == []: # if all vars in disp_rel are known
            is_valid = True if math.isclose(try_func,0,abs_tol=1e-7) else False
        else:
            print('adding missing parameter')
            solution = sym.solve(try_func, solve_for); # find value of missing parameter
            fk[var_solved[0]] = float(solution[0]); 
            is_valid = True; 
        return is_valid, fk

    def dispersion_rel(self, numeric=False):
        fk = self.wave_vector; 
        if numeric:
            kh2 = fk['kx']**2 + fk['ky']**2; bigK = kh2 + fk['kz']**2; 
            # Evaluate individual terms in dispersion relationship
            term_1 = self.medium['N2'] * kh2 / bigK; 
            term_2 = inertial_freq( self.medium['lat'])**2 * fk['kz']**2 / bigK;
            RHS = (term_1 + term_2)**(1/2); 
            relation = [[fk['omga'],RHS], ['omega','sqrt(N^2* kh^2/norm^2 + f^2 * kz^2 / norm^2)'] ];
        else:
            # Return ymbolic expression of dispersion relation
            relation = -self.lhs + self.rhs
        return relation

            

my_fk = {'omega':1e-4,'kx':0.000,'ky':2*np.pi/3e5,'kz':np.nan}
medium = {'lon':130, 'lat':20,'depth':200,'N2':4e-6}

niw1 = internal();
niw1.medium = medium; 
niw1.wave_vector = my_fk.copy(); # copy so my_fk won't change with class operations
print( niw1.wave_vector)

medium['N2'] = 0.9e-5; 
niw1.medium = medium; niw1.wave_vector = my_fk; 
print(niw1.wave_vector)
