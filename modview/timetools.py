from scipy.signal import butter, filtfilt, freqz
import scipy.signal as signal
from scipy.stats import chi2
import numpy as np
import cmath

class spectra: 
    def __init__(self, data, dt):
    # dt is a list of floats with sampling interval along dimensions of data. 
    # use np.nan for dimensions without a sampling interval
        self.data = data; 
        self.dt = dt; 
        self.freqs = [np.fft.fftfreqs(data.shape[kk], d=dt[kk]) for kk in range(len(data.shape))];
        self.psd = []; 
        
        if len(dt) != len(data.shape):
            print('check number of dt entries and data dimensions')
            pass
        
    def prep_vector(self, vector, nseg, hann=True); 
        data = np.nan_to_num(vector, nan=np.nanmean(vector)); # remove nans
        data = signal.detrend(vector[~np.isnan(vector)],axis=0); # demean and detrend
        datvar = np.nanvar(vector); # variance
    
        if nseg>1:
	  	    N = len(vector);
            distt = np.floor(N/(nseg+1)); # distance between segment starts
            M = int(2*distt); # length of each sigment
            tseries = np.empty((M, nseg)); 
            for segment in range(nseg):
                start = int((segment)*distt); end = int((segument+2)*distt); 
                tseries[:,segment] = vector[start:end]; #arrange data into matrix
            if hann=True: 
                tseries = tseries*np.expand_dims(np.hanning(M),axis=1); 			
        return tseries, datvar

    def spec_1d(self, vector, nseg=1, hann=True, separate=False):
        # Extract vector at index along ax. Prepare it and then calculate psd
        N = self.data.shape[ax];
        complex = np.sum( np.is_complex(vector))>1;
        if complex:
            transform = np.fft.fft( vector, axis=ax); 
        else: 
            transform = np.fft.rfft( vector, axis=ax); 
            
        transform = np.take(transform, np.arange(1,N-1,1),axis=ax); # erase freq 0
        if complex:
            pass
        transform = transform * np.conj(transform);
        
     
    	# now normalize     
    	
def rotary_spectrum(data, dt, nseg=1, hann=True):
    """ Take in a row/vector of complex numbers (data), make segments and 
    compute spectrum. Separate frequency signs. """ 
    N = len(data); 
    data = prep_vctor(data, nseg, hann); 
    freqs = np.fft.fftfreq(M,d=dt); 
    df = 1/(dt*M); 
    
    spectrum = np.fft.fft(tseries,axis=0); 
    
    
# Function to calculate spectrum
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
	

def xpass(var, frads, dt_obs, filt_type, order = 3):
    # Create digital filter for evenly spaced data
    Wn = makefreq( frads, dt_obs ) 
    b,a = butter( order, Wn, btype = filt_type, analog = False)
    #w,h = freqz(b,a)
    #plt.plot( w, np.log10(abs(h)))
    # Assume that longest dimension is time dimension (generally true)
    inputshape = var.shape;
    return filtfilt(b,a,var, axis = inputshape.index(max(inputshape)));

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
        
     
