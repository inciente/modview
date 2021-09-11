# Import necessary packages and classes 
import numpy as np
import xarray as xr # helps load and slice netcdf
import matplotlib.pyplot as plt
import scipy.signal as signal
from scipy.stats import chi2

# Load TAO SST from 8 N, 265 E as example 
datapath = 'data/sst8n95w_dy.cdf'
taodata = xr.open_dataset(datapath)

# Select subset of data with few gaps
taodata = taodata.sel(time=slice('2000-01-01','2008-12-31'));
# put data into 1D numpy array
taosst = np.squeeze( taodata['T_25'].values ); 

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
    data = np.nan_to_num( data, nan=np.nanmean(data)); # remove nans
    data = signal.detrend(data[~np.isnan(data)],axis=0); # remove mean and linear trend
    datvar = np.nanvar(data); 
    
    # Rearrange data into segments with 50% overlap
    if nseg!=1:
        distt = np.floor( N / (nseg+1)); # distance between segment starts
        M = int(2*distt); # length of each segment
        tseries = np.empty(( M , nseg));

        for segment in range(nseg):
            start = int((segment)*distt); end = int((segment+2)*distt); 
            tseries[:,segment] = data[start:end]; #arrange data into matrix

        if hann:
            tseries = tseries*np.expand_dims(np.hanning(M),axis=1); 
    else: 
        tseries = data;

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

# ---------- CALL SPECTRUM FUNCTION ----------
specsst = spectrum(taosst, 1, 6, True)

# ---------- PLOT RESULTS -----------
plt.loglog( specsst['freqs'], specsst['spectrum'], label='spectrum')
plt.semilogy( [3e-1, 3e-1], specsst['errbar']/20,linewidth=3, label='errorbar')
plt.grid(True,which='major'); plt.legend()
plt.ylabel('Spectral density [C$^2$ / cpd]'); plt.xlabel('Frequency [cpd]'); 
plt.show()
