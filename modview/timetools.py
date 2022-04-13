from scipy.signal import butter, filtfilt, freqz
import scipy.signal as signal
from scipy.stats import chi2
import numpy as np
import cmath, datetime
import pandas as pd
from abc import ABC, abstractmethod
import xarray as xr

def prep_data(data_mat, axis=0, detrend=False):
    avg_real = np.nanmean(np.real(data_mat)); 
    avg_imag = np.nanmean(np.imag(data_mat)); 
    data_mat = np.nan_to_num( np.real(data_mat), nan=avg_real) \
                  + 1j*np.nan_to_num( np.imag(data_mat), nan=avg_imag);
    if detrend: 
        data_mat = signal.detrend( np.real(data_mat), axis=axis) \
                  + 1j*signal.detrend( np.imag(data_mat), axis=axis); 
                  # output will always be imaginary
    return data_mat
    
def nan_to_mean(vector):
    checknan = np.isnan(vector); 
    cavg = np.nanmean(vector); 
    vector[checknan] = cavg;
#    nuvector = np.nan_to_num(float(np.real(vector)), nan=float(np.real(cavg))); 
    return vector

def segment_vector(vector, nseg, hann=False, detrend=False, axis=0):
    ''' Take a vector (for example a time serise) and separate it into nseg segments with a 50% overlap.
    Options allow to add detrend the entire vector and/or scale windows using hanning window.
    Axis is the index along which data will be segmented. Other axes are taken to be 
    separate observations and same operation is repeated for those.'''
    # Replace nans for mean of vector 
    data = prep_data( vector, axis, detrend ); 
    datvar = np.nanvar(data,axis=axis) # compute data variance
    if nseg>1:
        # SEPARATE DATA INTO CHUNKS
        N = len(vector); # length of input
        distt = np.floor(N/(nseg+1))
        M = int(2*distt); # length of each resulting segment 
        tseries = np.empty((M,nseg)); 
        for segment in range(nseg): # save chunks in tseries matrix 
            start = int((segment)*distt); end = int((segment+2)*distt); 
            tseries[:,segment] = data[start:end]; #arrange data into matrix
        if hann==True: 
            tseries = tseries*np.expand_dims(np.hanning(M),axis=1);   			
        return tseries, datvar # return segmented and windowed data
    else:
        return data, datvar # return single segment of (detrended) data

def power_from_mat(data_mat, axis=0, detrend=False):
    # Enter a 2D array and return the spectra of columns (ax=0) or rows (ax=1)
    N = data_mat.shape[axis]; # number of elements for fft
    M = data_mat.shape[np.abs(axis-1)]; # number of realizations
    data_mat = prep_data(data_mat, axis, detrend); # remove nans, detrend

    comp = np.sum( np.iscomplex(data_mat))>1; # check if complex
    if comp:    
        transform = np.fft.fft( data_mat, axis=axis); 
    else: 
        transform = np.fft.rfft( np.real(data_mat), axis=axis); 
    # make sure line below works for both fft and rfft 
    #transform = np.delete( transform, 0, axis=axis); # remove freq=0 
    transform = np.real(transform * np.conj(transform)); # calc power       
    return transform

def spectrum_1D(arr, dt, axis, nseg=1, hann=False, **kwargs):
    arr = np.squeeze(arr); # avoid shape issues
    N = arr.shape[axis]; df = 1/(N*dt); 
    data, datvar = segment_vector(arr, nseg, hann, axis=axis);
    datvar = np.expand_dims(datvar, axis=axis);  
    power = power_from_mat(data, axis=axis); 
        
    #if len(power.shape)>1 and power.shape[1]>1:
    #    power = np.mean( power, axis=1); 
    powvar = np.expand_dims( np.sum(power,axis=axis)*df, axis=axis); # integrate spectra 
    power = power / powvar * datvar; # normalize
    return power

def spectral_error(degsfred, cl=0.95):
    err_low = chi2.isf(cl+(1-cl)/2, degsfred); 
    err_hi = chi2.isf((1-cl)/2, degsfred); 
    err_bar = np.array([err_low, err_hi]); 
    return err_bar          


class spectra: 
    def __init__(self, data, dt):
    # dt is a list of floats with sampling interval along dimensions of data. 
    # use np.nan for dimensions without a sampling interval
        self.data = data; 
        self.dt = dt;
        self.psd = dict(); 
        if len(dt) != len(data.shape):
            print('check number of dt entries and data dimensions')
            pass
        
    def prep_vector(self, vector, nseg, hann=True, ax=0):
        data = np.nan_to_num(vector, nan=np.nanmean(vector)); # remove nans
        data = signal.detrend(data,axis=ax); # demean and detrend
        datvar = np.nanvar(vector); # variance
    
        if nseg>1:
            N = len(vector);
            distt = np.floor(N/(nseg+1)); # distance between segment starts
            M = int(2*distt); # length of each sigment
            tseries = np.empty((M, nseg)); 
            for segment in range(nseg):
                print(segment)
                start = int((segment)*distt); end = int((segment+2)*distt); 
                tseries[:,segment] = data[start:end]; #arrange data into matrix
            if hann==True: 
                tseries = tseries*np.expand_dims(np.hanning(M),axis=1);
        else: 
            tseries = data;    			
        return tseries, datvar

    def get_power(self, data_mat, ax=0):
        N = self.data.shape[ax];
        comp = np.sum( np.iscomplex(data_mat))>1;
        if comp:
            transform = np.fft.fft( data_mat, axis=ax); 
        else: 
            transform = np.fft.rfft( data_mat, axis=ax); 
        transform = np.delete( transform, 0, axis=ax); # remove freq=0
        transform = np.real(transform * np.conj(transform)); # calc power
        return transform
        
    def spec_vec(self, vector, dt, trans_ax, nseg=1, hann=True, freqs=False):
        N = len(vector); df = 1/(N*dt); 
        data, datvar = self.prep_vector(vector, nseg, hann, ax=trans_ax); 
        power = self.get_power(data, ax=trans_ax); 
        print(power.shape)
        if len(power.shape)>1 and power.shape[1]>1:
            power = np.mean( power, axis=1); 
        powvar = np.sum(power,axis=trans_ax)*df; 
        power = power / powvar * datvar; # normalize
        return power
    
    def err_bars(self, degfred):
        err_low = chi2.isf(0.975, degsfred); 
        err_hi = chi2.isf(0.025, degsfred); 
        err_bar = np.array([err_low, err_high]); 
        return err_bar

    def calc_spec(self, ax, indices, nseg=1, hann=True):
        pass
    
    
# Function to calculate spectrum
'''
def spectrum(data, dt, nseg=1, hann=True):
    """
    Separate the timeseries "data" into "nseg" separate segments, calculate their
    average frequency spectra, and return a dictionary with all variables necessary 
    for analysis. 
    - - "dt" is the time interval between datapoints (arbitrary units).
    - - "data" must be a numpy array with shape (N,1).
    - - Output has frequency units cycles/[units of dt].
    """
    N = len(data); 
    tseries, datvar = prep_vector(data, nseg, hann); 
    # Frequency info before getting spectrum
    freqs = np.fft.rfftfreq( M, d=dt); 
    freqs = freqs[1:]; # remove freq = 0 (mean)
    df = 1/(dt*M);
    
    spectrum = np.fft.rfft(tseries, axis=0); # real fft
    spectrum = spectrum*np.conj(spectrum); 
    spectrum = np.real( np.delete( spectrum, 0, axis=0) ); # remove row for freq=0
    
    avg_spectrum = np.mean( spectrum, axis=1); # average across segments
    varinspec = np.expand_dims( np.nansum( avg_spectrum*df, axis=0),axis=0 ); #area under spectrum
    avg_spectrum = avg_spectrum / varinspec * datvar; # normalized 
    
    degsfree = 36*nseg/19; # degrees of freedom of avg_spectrum
    err_low = chi2.isf( 0.975, degsfree);
    err_high = chi2.isf(0.025, degsfree); 
    err_bar = np.array( [err_low, err_high])
    
    result = {'spectrum':avg_spectrum,'freqs':freqs,'df':df,'errbar':err_bar};
    return result
'''	

def powerspec(timeseries, dt, timeax=0, hann=True):
	"""
	Calculate the power density spectrum of object "timeseries" (must be real)
	"""

	# Store variance of each segment:
	normvar = np.expand_dims( np.var( detrended, axis=timeax), axis=timeax); 
	df = 2*np.pi/(detrended.shape[timeax]*dt); 
	
	if hann:
		detrended = detrended*np.hanning(detrended.shape[timeax]); # hanning window
	
	spectrum = np.fft.rfft(detrended, axis=timeax); #real fft
	spectrum = spectrum*np.conj(spectrum); #squared
	spectrum = np.delete( spectrum, 0, axis=timeax); # delete row for freq=0
	
	varinspec = np.expand_dims( np.nansum( spectrum*df, axis=timeax), axis=timeax); # area under spectrum
	spectrum = (spectrum/varinspec)*normvar; # normalize for variance
	
	return np.real(spectrum)
	

def xpass(var, frads, dt_obs, filt_type, order = 3, axis=None):
    # Create digital filter for evenly spaced data
    Wn = makefreq( frads, dt_obs ) 
    b,a = butter( order, Wn, btype = filt_type, analog = False)
    #w,h = freqz(b,a)
    if axis is None:
        inputshape = var.shape;
        axis = inputshape.index(max(inputshape));
    
    #plt.plot( w, np.log10(abs(h)))
    # Assume that longest dimension is time dimension (generally true)
    inputshape = var.shape;
    return filtfilt(b,a,var, axis=axis);

def makefreq(frads, dt_obs):
    # Create frequency input for xpass based on:
    # frads  - desired frequency in radians / second
    # dt_obs - sampling rate in datapoints / second
    
    # Scipy.signal tools take frequencies in units (half cycles / sample). 
    normT = 1/dt_obs; # seconds per sample
    Wn = frads / np.pi ; # convert to half-cycles / second. 
    Wn = Wn * normT ; # half-cycles per sample 
    return Wn     
    
def nearf(lat,width = [0.8,1.2]):
    # Find band-limit frequencies for the near-inertial band
    f = 4*np.pi*np.sin(lat/180*np.pi)/24/3600; # inertial freq in rads/second
    return f*np.array(width)


def get_timenum(var_xr, units='seconds'):
    timenum = pd.to_datetime(var_xr.time.values);
    timenum = (timenum - timenum[0]).total_seconds();
    timenum = np.expand_dims( timenum.to_numpy(), axis=1); 
    if units in ['hours','days']:
        timenum = timenum/3600;
    if units == 'days':
        timenum = timenum / 24;
    return timenum


class Hebert_1994:
    ''' Backrotation following the appendix of Hebert and Moum (JPO, 1994).
    Mostly meant to be used with shipboard adcp data, but can be used in other
    settings too. '''
    def __init__(self, var_xr, seg_dt, freqs=np.linspace(0.95,1.05,4)):
        var_xr['timenum'] = ('time', np.squeeze(get_timenum(var_xr)) ); # unique time vector
        self.ship_dat = var_xr;
        self.seg_dt = seg_dt; # segment duration in seconds
        self.freqs = freqs; # to be multiplied by local f for each segment
        self.results = dict();
        print('about to get amplitudes')
        self.get_amplitudes()

    def backrotate(self, dataset, variable, frads, forward=False, amp=None):
        timephase = frads * dataset['timenum']; # use universal time vector
        if isinstance(variable, str):
            # prepare dataset[variable] for backrotation
            variable = dataset[variable]; 
            variable = variable - variable.mean(dim='time', skipna=True);

        if forward: 
            if isinstance(amp, xr.DataArray):
                amp = amp - amp.mean(dim='time', skipna=True);
            else: 
                amp = np.expand_dims(amp,axis=0);
                timephase = timephase = np.expand_dims(timephase,axis=1); 
            backrotated = amp * (np.cos(timephase) + 1j*np.sin(-timephase));
        else:
            variable = variable - variable.mean(dim='time',skipna=True);
            backrotated = variable * ( np.cos(timephase) + 1j*np.sin(timephase) );
        return backrotated

    def variance_fraction(self, A,B,freqs, segment_data):
        # Calculate variance in data for a list of amplitudes and frequencies
        u_dat = segment_data['u'] - segment_data['u'].mean(dim='time',skipna=True);
        v_dat = segment_data['v'] - segment_data['v'].mean(dim='time',skipna=True); 
        denom = u_dat.var(skipna=True) + v_dat.var(skipna=True);  
        fractions = [] # save data here
        for kk in range(len(freqs)):
            amp_u = A[kk]; amp_v = B[kk]; frads = freqs[kk]; 
            u_rec = self.backrotate(segment_data,'u',frads,forward=True,amp=amp_u);
            v_rec = self.backrotate(segment_data,'v',frads,forward=True,amp=amp_v);
            # Calculate MSE
            numer = np.abs( u_dat - u_rec )**2 + np.abs( v_dat - v_rec )**2;
            numer = np.nanmean(numer); 
            fractions.append(numer);
        fractions = 1 - np.array(fractions) / numer; 
        return fractions

    def get_amplitudes(self):
        # Run through segments and calculate amplitudes
        Ti = pd.to_datetime( self.ship_dat.time.values[0] ); 
        Tf = Ti + datetime.timedelta(seconds=self.seg_dt); # segment limits
        A = []; B = []; omega = []; Times = []; # to save answers
        var = []; local_f = []; 
        while Tf < self.ship_dat.time.values[-1]:
            segment = self.ship_dat.sel(time=slice(Ti,Tf)); # apply time limits
            f_test = nearf(segment['latitude'].mean().values, self.freqs); # ni-range
            local_f.append( nearf(segment['latitude'].mean().values, 1));
            # ----- save fit information for many frequencies
            amps_u = []; amps_v = [];
            for kk in range(len(f_test)):
                rotated_u = self.backrotate(segment,'u', f_test[kk]); 
                rotated_v = self.backrotate(segment,'v', f_test[kk]); 
                # save complex amplitudes)
                u_amp = 2*np.nanmean( np.abs(rotated_u.values),axis=0); 
                u_phase = 2*np.nanmean(rotated_u.values,axis=0)/u_amp
                amps_u.append( u_amp*(np.cos(u_phase) + 1j*np.sin(u_phase) )); 
                # now for v
                v_amp = 2*np.nanmean( np.abs(rotated_v.values),axis=0);
                v_phase = 2*np.nanmean( rotated_v.values,axis=0)/v_amp;
                amps_v.append( u_amp*(np.cos(v_phase) + 1j*np.sin(v_phase)  ))
            # ----- compute fraction of variance in each fit
            var_fracs = self.variance_fraction( amps_u, amps_v, f_test, segment);
            # keep best fit only
            best_index = np.argmax(var_fracs);
            A.append(amps_u[best_index]); B.append(amps_v[best_index]);
            omega.append(f_test[best_index]); 
            var.append( var_fracs[best_index] ); 
            Ti += datetime.timedelta(seconds=self.seg_dt/2); # for next segment
            Tf += datetime.timedelta(seconds=self.seg_dt/2);# 50% overlap
            Times.append(Ti); # save middle point in segment
        
        self.results = {'A':A,'B':B,'omega':omega,'Tm':Times,'f':local_f}
    
    def remake_data(self,comp='u', dt=7200):
        results = self.results
        if comp == 'u':
            amp = results['A'];
        elif comp == 'v':
            amp = results['B'];
        recreate = xr.Dataset( coords={'pressure':self.ship_dat['depth'].values,
                        'time':np.squeeze(results['Tm']) } );
        recreate['amp'] = xr.DataArray( data=np.asarray(amp),
                dims=('time','pressure'))
        recreate['omega'] = xr.DataArray( data=results['omega'], dims='time')
        recreate = recreate.interp(time=self.ship_dat.time); 
        # Alternative way of computing this shit
        timephase = self.ship_dat['timenum']*recreate['omega']; 
        recreate['timephase'] = timephase
        u_rec = recreate['amp']*( np.cos(timephase) + 1j*np.sin(\
                timephase));
        
       # ------- time vector to interpolate onto        
       # tot_duration = (results['Tm'][-1] - results['Tm'][0]).total_seconds();
       # t_off = np.arange(0, tot_duration, dt);
       # new_time = [results['Tm'][0] + datetime.timedelta(seconds=kk) \
       #         for kk in t_off]
       # recreate = recreate.interp(time=new_time); 
       # timephase = get_timenum(recreate).transpose()*recreate['omega'].values; 
       # timephase = timephase.transpose();
       # recreate['timephase'] = (('time'), np.squeeze(timephase)); 
       # # apply rotation
       # u_rec = recreate['amp'].values*( np.cos( timephase ) + \
       #         1j*np.sin( timephase ) );
        u_rec = xr.DataArray( data=u_rec.transpose(), 
                dims=('pressure','time'), coords={\
                'pressure':self.ship_dat['depth'].values,'time':self.ship_dat.time})
        return u_rec, recreate

    def NIKE(self,z0,z1,dt=3600*12):
        u_fake = self.remake_data('u',dt=dt);
        v_fake = self.remake_data('v',dt=dt); 

    def u(self, v=False):
        if v:
            u = self.ship_dat['v'].transpose()
        else:
            u = self.ship_dat['u'].transpose();
        return u

def CompDemod(variable, tvec, frads, dt_obs):
    # Data in variable taken to be periodic at freq + other things.
    timevec = np.array([kk for kk in range(len(tvec))]);
    timevec = np.expand_dims(timevec, axis =1 );
    # make sin wave with the right frequency
    period = 2*np.pi/frads; points_in_period = period*dt_obs; 
    timevec = -timevec/points_in_period*2*np.pi
    #timevec = -2*freq*timevec; # turn time into phase
    Yvar = variable*( np.cos(timevec) + 1j*np.sin(timevec+np.pi/2));

    # Get rid of nans before using filter
    Yvar = np.nan_to_num(Yvar); 
    Yvar = xpass( Yvar, 0.75*frads, dt_obs, 'low',order = 4); # filter
    # Mask data originally unavailable
    #Yvar = Yvar.transpose(); 
    Yvar[ np.isnan( variable )] = np.nan; 

    Amp = 2*np.abs(Yvar); cycphase = 2*Yvar/Amp; 
    return Amp

     
