import numpy as np; import xarray as xr; import pandas as pd; 
import sympy as sym; import gsw; import mapper; import datetime
import scipy.interpolate as scinterp
from sympy.utilities.lambdify import lambdify, implemented_function
import importlib
importlib.reload(mapper)

def lon_to_x(lon, lat):
    x = lon*110e3*np.cos(np.mean(lat)); x = x - np.mean(x);
    return x

def lat_to_y(lat): 
    y = lat*110e3; y = y - np.mean(y);
    return y
    
def inertial_freq(lat):
    f = 4*np.pi*np.sin(lat/180*np.pi)/24/3600;
    return f
    
def hydrostasy(rho, zlevels, ssh=0, z_req=np.nan, p_axis=0, p_ref = 0):
    # Take in n-dimensional object rho and integrate pressure 
    # ssh represents local values at locations of density data
    # currently only works for 1d objects of density
    dz = np.diff(zlevels)
    pressure = np.empty(rho.shape); 
    pressure[:] = np.nan;
    pressure[0] = (ssh + zlevels[0])*g*min(rho); 
    
    pressure[1:] = [dz[i]*9.81*0.5*(rho[i-1]+rho[i]) for i in range(1,len(pressure))];
    return pressure

def wkb_stretch(z, N2, MLD):
    belowML = z>MLD; 
    N = np.sqrt(N2); 
    sratio = np.nanmean(N[belowML])/N; 
    sratio = pd.DataFrame(sratio); # to get rolling mean
    sratio = sratio.rolling(10, min_periods=1).mean().to_numpy(); 
    
    vel_factor = np.sqrt(sratio); # velocity scaling
    
    tot_depth = z[~belowML]; # stretch vertical coordinate
    stretched = np.expand_dims( np.append(0, np.diff(z)), axis=1); 
    stretched = np.cumsum( stretched[belowML]/sratio[belowML])+MLD; # vertical integration
    stretched = np.append(tot_depth, stretched); 
    return vel_factor, stretched
   

# WKB STRETCHING AND SCALING FACTORS
argo_path = '/mnt/sda1/SciData/ARGO/RG_ArgoClim_33pfit_2019_annual.nc';
def argo_wkb(lon,lat,adcp_dat,zcoord='depth'):
    argo = xr.open_dataset(argo_path);
    argo = argo.sel(LONGITUDE=134, method='nearest');
    argo = argo.sel(LATITUDE=14, method='nearest');
    # Save everything necessary to perform wkb stretching
    ctd_prof = {'T':argo['ARGO_TEMPERATURE_MEAN'],
                'S':argo['ARGO_SALINITY_MEAN'],
                'p':argo['PRESSURE']}
    ctd_prof['CT'] = gsw.CT_from_t(ctd_prof['S'],ctd_prof['T'], ctd_prof['p'])
    ctd_prof['N2'], ctd_prof['p_mid'] = gsw.Nsquared(ctd_prof['S'],ctd_prof['CT'],ctd_prof['p']);
    # Now define wkb properties
    wkb_prof = {'z':adcp_dat[zcoord]};
    # n2 at adcp bins:
    wkb_prof['N2'] = np.interp(wkb_prof['z'],ctd_prof['p_mid'], ctd_prof['N2'])
    wkb_prof['u_scale'], wkb_prof['z_stretch'] = wkb_stretch( \
                wkb_prof['z'],wkb_prof['N2'], 30 );
    return wkb_prof, ctd_prof

def adcp_to_wkb(adcp_block, wkb_prof, dz=10, axis=0):
    ''' Take in xr.DataArray representing ADCP measurements and 
    1) Scale velocities using wkb_prof['u_scale']
    2) Interpolate data from stretched coordinate wkb_prof['z_stretch']
     ---  onto constantly-spaced version
    '''
    block = adcp_block * wkb_prof['u_scale'];
    z_wkb = wkb_prof['z_stretch'];
    z_wkb_c = np.arange(min(z_wkb), max(z_wkb)+1, dz);
    
    # Now use inteprolator to go from z_stretch to z_constant_dz
    interpolator = scinterp.interp1d( wkb_prof['z_stretch'], block.values,
            kind='slinear', axis=axis);
    block = xr.DataArray(data=interpolator(z_wkb_c),
            dims=['depth','time'],
            coords={'depth_wkb':('depth',z_wkb_c), 'time':('time',block.time) } );
    return block

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
        
def ray_tracing(wave, medium, timevec):
    i0 = wave.fk; # initial properties of wave
    p0 = wave.position; # initial wave location
    drdt = wave.group_vel(); # can method be saved?
    # need to invoke c_g functions and run them through a generic solver
    dkdt
    return 
             
class internal_wave:
    def __init__(self, lat, t0=0,N2=1e-5):
        self._fk = dict();  # vector [freq, kx, ky, kz].
        self._Nsquared = N2;
        self.position = {'lat':lat};
        self.t0 = t0; # when was it generated? 
        self.lat = lat;
        # Add new method that extracts N2 from a stratification of geostrophic object. 
    
    @property
    def Nsquared(self):
        return self._Nsquared
    
    @Nsquared.setter
    def Nsquared(self, geost_field):
        # Set value to argument n2
        # Access methods and properties of geostrophic field to extract N2 at (lon, lat).
        n2 = mapper.var_around( geost_field.N2, self.position); 
        # is position in the right format?
        self._Nsquared = n2;
        
    @property
    def fk(self):
        return self._fk
    
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
            if try_func == 0:
                self._fk = wave_properties;
            else:
                print('wave properties not compliant with dispersion relationship')
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
    
    def fk_sym(self):
        wave_properties = [(sym.symbols(key),value) for key, value in self._fk.items()]
        return wave_properties
    
    def fk_vals(self): 
        wave_properties = self.fk;
        kx = wave_properties['kx']; ky = wave_properties['ky']; 
        kz = wave_properties['kz']; omega = wave_properties['omega'];
        return omega, kx, ky, kz

    def set_position(self, lon, lat, depth):
        self.position = [lon, lat, depth]; 
        self.lat = lat; 
        
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
        if kunze==False:
            disp_rel = self.dispersion_relation();
        else:
            disp_rel = self.kunze_omega_0('none',numeric=False); 
        
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
            print(ratio)
            return False
        
    def kunze_omega_0(self, geost_field, which_terms=[1,2,3], numeric=True):    
        # input necessary: dU/dz, dV/dz, dU/dy, dV/dx
        if self.is_niw() == False:
            print('wave is not near-inertial')
            return
        else:
            f_loc = inertial_freq(self.lat); 
            if numeric:
                dudz, dvdz = geost_field.thermal_wind(self.position, point=True); 
                # write method to get vorticity from density and ssh (their laplacian)
                dvdx, dudy = geost_field.vorticity(self.position, point=True); 
                f_eff = f_loc + 0.5*(dvdx-dudy);
                omega, kx, ky, kz = self.fk_vals(); 
                kh_squared = kx**2 + ky**2; 
                # make list with terms in omega_0
                om_terms = [0, f_eff, # terms 0, i
                       self.Nsquared*kh_squared/(2*f_loc*kz**2), # term ii
                       1/kz*( dudz*ky - dvdz*kx), # term iii
                       1j*( vort/2/f_loc/kz*(dudz*kx + dvdz*ky)),  # term iv
                       self.Nsquared/(2*f_loc**2*kz**2)*(dudx*ky**2\
                       - (dudy + dvdx)*kx*ky + dvdy*kx**2 )]; # term v   
                # intrinsic frequency
                omega_0 = np.sum([om_terms[kk] for kk in range(6) if kk in which_terms]);
                # account for doppler shift
                #omega = om_terms[1] + om_terms[2] + om_terms[3] \
                #        + kx*u + ky*v; # equation 14, assume w = 0
            else:
				# use sympy
                wave_p = self.fk_sym(); 
                zeta, dudz, dvdz, N2 = sym.symbols('zeta dudz dvdz N2');
                kx = wave_p[1][0]; ky = wave_p[2][0]; kz = wave_p[3][0]; 
                om_terms = [0, f_loc + zeta/2, # term i
                           N2*(kx**2 + ky**2)/(2*f_loc*kz**2), # term ii
                           1/kz*(dudz*ky - dvdz*kx)]; # term iii
                omega_0 = omega_0[0];
                for kk in [1,2,3]:
                    if kk in which_terms:
                        omega_0 += om_terms[kk];

            return omega_0
                 
    def solve_kunze(self, geostr, timevec):
		# use wave properties as initial condition and save r, k 
		# solutions at times in timevec 
		# invoke setup_kunze
        k_0 = self.fk_vals(); 
        r_0 = self.position; # lon, lat, depth
		 
		# use scipy ode solvers to integrate dk/dt        
    
    def internal_tide(self,component=['M2','D1']): 
        # Set f_k to represent an internal tide
        pass
        
    def plane_wave(self, as_dates=True, offset=0, **kwargs):
        # Use wave characteristics to create a multi-dimensional sinusoidal wave. 
        # Read kwargs as if it were a coords dictionary and infer desired dimensions 
        pos_dims = ['time','x','y','z']; 
        if 'time' in kwargs:
            timevec = kwargs['time'].copy(); 
            if as_dates:
                num_time = [(timevec[kk] - timevec[0]).total_seconds() \
                           for kk in range(len(timevec))];
                kwargs['time'] = num_time;   
                
        omega, kx, ky, kz = self.fk_vals(); # wave properties
        sine_vals = [];  
        for key, value in kwargs.items():
            print(key)
            if key == 'x': # identify dimension 
                phase_num = kx;
            elif key == 'y':
                phase_num = ky;  
            elif key == 'z':    
                phase_num = kz; 
            elif key == 'time':
                phase_num = omega;
            dim_phase = np.array(value)*phase_num 
            sine_vals.append( np.cos(dim_phase) + 1j*np.sin(dim_phase) );
        # Manipulate dimensions for multiplication
        ndims = len(sine_vals); 
        dim_lens = []; # desired shape of out
        for kk in range(ndims):
            dim_lens.append(len(sine_vals[kk])); 
            expansion = tuple([jj for jj in range(ndims) if jj!=kk]) # remove dimension of data
            sine_vals[kk] = np.expand_dims(sine_vals[kk], expansion); 
        # Calculate wave function
        wave_vals = np.prod(sine_vals);  
        return wave_vals
        
class geostrophic_field:
    def __init__(self, CT, SAL, P, dims='both', ssh='none'):
        # all inputs are xr.dataarrays. SSH may have different coordinates
        self.CT = CT; 
        self.SAL = SAL; 
        self.P = P;
        self.density = gsw.density.rho(self.SAL, self.CT, self.P); 
        self.p_axis = self.find_axis(self.SAL, self.P); 
        self.psi = self.streamfunction();
        self.ssh = ssh;
    
    def streamfunction(self): 
        p_ax = self.p_axis;        
        p_mat = self.upsize_p(self.CT);
        DH = gsw.geo_strf_dyn_height(self.SAL.values, self.CT.values, p_mat, axis=p_ax);
        DH = xr.DataArray(data=DH, dims=self.CT.dims, coords=self.CT.coords); 
        return DH
    
    def var_around(self, xr_obj, location, edge=1): 
		# Take in coordinates, and find u,v at that spot and neighboring points.
		# location is [lon, lat, depth]
        subsel = mapper.block_around(xr_obj, location, edge=edge); 
        return subsel
        
    def lat_lon_deriv(self,variable, location):    
        var_x = xr.DataArray(data=np.nan, coords=location); 
        var_y = xr.DataArray(data=np.nan, coords=location);
        if 'lat' in variable.dims:
            var_y = variable.differentiate(coord='lat')/110e3;
        if 'lon' in variable.dims:
            var_x = variable.differentiate(coord='lon')/110e3 \
                         / np.cos(np.radians(location['lat'])); 
        return var_x, var_y
        
    def eta_to_u(self, location, edge=1, use='psi', point=False):
		# use is set to psi (streamfunction) or 'ssh'. 
		# in both cases, we take horizontal derivatives to calculate u,v
        c0 = 9.81 / inertial_freq(location['lat']); 
        if use=='ssh': 
            del location['p'];
            myvar = self.ssh;  
        elif use=='psi': # use appropriate function
            myvar = self.psi; 
        elif use=='rho':
            myvar = self.density; 
            c0 = -c0 / 1023; 
        if location=='none':
            var_here = myvar; 
        else:
            var_here = self.var_around( myvar, location, edge=edge); 
        # Differentiate over x and y.     
        [v_here, u_here] = self.lat_lon_deriv(var_here, location); 
        u_here = -c0 * u_here; 
        v_here = c0 * v_here; 
        if point: # evaluate at single point?
            u_here = mapper.middlepoint(u_here.values); 
            v_here = mapper.middlepoint(v_here.values);          
        return u_here, v_here        
              
    def flow_at(self, location, edge=1, point=False):
        # Calculate u,v from stratification
        u_strat, v_strat = self.eta_to_u(location, edge, use='psi', point=point); 
        if self.ssh != 'none':
            # Calculate barotropic component
            u_bar, v_bar = self.eta_to_u(location, edge, use='ssh'); 
            u_bar = mapper.match_latlon(u_bar, u_strat); 
            v_bar = mapper.match_latlon(v_bar, v_strat); 
        else:
            u_bar = 0; v_bar = 0;
        u = u_bar + u_strat; 
        v = v_bar + v_strat; 
        if point==True:
            u = mapper.middlepoint(u);
            v = mapper.middlepoint(v); 
        return u, v
		       
    def thermal_wind(self, location, edge=1, point=False):
        dudz, dvdz = self.eta_to_u(location, edge, use='rho');
        if point:
            dudz = mapper.middplepoint(dudz.values); 
            dvdz = mapper.middlepoint(dvdz.values); 
         
        return dudz, dvdz
    
    def vorticity(self, location, edge=[0,2,2], point=False):
        u,v = self.flow_at(location,edge); 
        dvdx = v.differentiate(coord='lon')/110e3/np.cos(np.radians(location['lat']));
        dudy = u.differentiate(coord='lat')/110e3;  
        if point: # return value in the middle of block
            dvdx = mapper.middlepoint(dvdx.values); 
            dudy = mapper.middlepoint(dudy.values); 
        vort = dvdx - dudy
        return vort

         
    def pressure(self, use_ssh=False):
        if use_ssh: # barotropic component
            eta = self.ssh.sel(lat=self.density.lat, method='nearest'); 
            eta = eta.sel(lon=self.density.lon, method='nearest'); 
        else:
            pass
        dens_pres = 9.81*self.density.integrate(coord='p'); # effect of stratification
        pres = 9.81*1023*eta + dens_pres; 
        return pres
      
    def upsize_p(self, var):
        # Repeat pressure vector to match the shape of var
        p_axis = self.find_axis(var, self.P);
        new_p = np.ones(var.shape); 
        # Reshape pressure to populate columns in new_p volume
        ax_to_create = [kk for kk in range(len(var.shape)) if kk != p_axis]; 
        p_col = np.expand_dims( self.P.values, tuple(ax_to_create) );
        new_p = new_p * p_col; 
        return new_p
        
    def subsel(self, variable, location, loc_format='coords'):  
        # location  is a dictionary of lists
        if loc_format == 'coords':
            var_subsel = variable; # save XR object to be sliced
            for key in location:
                var_subsel = var_subsel.sel(key=location[key], method='nearest');   
    def find_axis(self, var, vec):
        check = [var.shape[kk] == len(vec) for kk in range(len(var.shape))]; 
        axis = [i for i, x in enumerate(check) if x]; 
        if len(axis)>1:
            print('dimensions agree along more than one axis')
            return
        else:
            return axis[0]

    def currents(self, location, loc_format='coords', calc_gradients=False):
        # location must be a dictionary with format {'lat':[12.3,12.5],'lon':[145,148],
        # 'p':[250, 300]; 

        # Now that we have extracted only the necessary portion of data, we need to
        # get spatial derivatives of the stream function to return uv. 
        U = -sf_subsel.differentiate(coord=lat)/110e3;
        V = sf_subsel.differentiate(coord=lon)/110e3/ \
               np.cos(np.nanmean(location['lat'])/180*np.pi); 
        return U, V

def kunze_1985(niws, geost, slow_time, fast_time):
    # Add functions necessary to solve dkdt and drdt using the properties of 
    # internal waves. 
    # Add two iterative cycles: one for slow time (update geostrophic fields), 
    # another for fast time (NIW solver)
    
    # Import terms of geostrophic shear: symbols and numerical grids.
    u, v = geost.flow_at(location='none', point=False); 
    dvdx, dudy = geost.vorticity(location='none', point=False);
    dudz, dvdz = geos.thermal_wind(location='none', point=False); 

    # Consider allowing to use sets of niws, so computing dedicated to geos can be
    # reused for many niws. 
    # niws object will have to include information about the energy of internal waves and when they should be seeded.
    
        
    # Having loaded the entire geostrophic field, we can now use mapper.var_around for each wave
    
     
    
    # Allow access to a method that calculates omega according to kunze
    # For each wave
    # ------ calculate omega (starting condition)
    # ------ dkdt : evaluate grad(omega) at niw.fk (each timestep)
    # ------ drdt : evaluate group velocity (each timestep) 
    # ------ solver: allow access to niw.fk
    




