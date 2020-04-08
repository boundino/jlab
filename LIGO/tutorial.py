#----------------------------------------------------------------
# Load libraries
#----------------------------------------------------------------

import numpy as np
import math
from scipy import signal
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, iirdesign, zpk2tf, freqz

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import h5py

# LIGO-specific readligo.py 
import readligo as rl

#----------------------------------------------------------------
# Load LIGO data
#----------------------------------------------------------------

fn = 'data/H-H1_LOSC_4_V1-1126259446-32.hdf5'
strain, sampletimes, chan_dict = rl.loaddata(fn, 'H1')

#----------------------------------------------------------------
# Show LIGO strain vs. time
#----------------------------------------------------------------

plt.figure()
plt.plot(sampletimes, strain)
plt.xlabel('Time (s)')
plt.ylabel('strain')
plt.title('Advanced LIGO strain data near GW150914')
plt.savefig('plots/GW150914_strain.png')

#----------------------------------------------------------------
# Obtain the power spectrum density PSD / ASD
#----------------------------------------------------------------

psdx, psdfreqs = mlab.psd(strain, Fs = 4096, NFFT = 4096)
psd = interp1d(psdfreqs, psdx)

plt.figure()
plt.loglog(psdfreqs, np.sqrt(psdx),'r',label='H1 strain')
plt.axis([10, 1500, 1e-24, 1e-19])
plt.ylabel('ASD (strain/rtHz)')
plt.xlabel('Freq (Hz)')
plt.title('Advanced LIGO strain data near GW150914')
plt.savefig('plots/GW150914_ASDs.png')

#----------------------------------------------------------------
# Whitening data
#----------------------------------------------------------------

def whiten(strain, interp_psd, dt):
    Nt = len(strain)
    freqs = np.fft.rfftfreq(Nt, dt)
    hf = np.fft.rfft(strain)
    white_hf = hf / (np.sqrt(interp_psd(freqs) /dt/2.))
    white_ht = np.fft.irfft(white_hf, n=Nt)
    return white_ht

white_data = whiten(strain, psd, sampletimes[1]-sampletimes[0])

plt.figure()
plt.plot(sampletimes, white_data)
plt.xlabel('Time (s)')
plt.ylabel('strain')
plt.title('Advanced LIGO strain data near GW150914 with whitening')
plt.savefig('plots/GW150914_strain_whiten.png')

#----------------------------------------------------------------
# Bandpass filtering
#----------------------------------------------------------------

bb, ab = butter(4, [20.*2./4096, 300.*2./4096], btype='band')
white_data_bp = filtfilt(bb, ab, white_data)

plt.figure()
plt.plot(sampletimes, white_data_bp)
plt.xlabel('Time (s)')
plt.ylabel('strain')
plt.title('Advanced LIGO strain data near GW150914 filtered')
plt.savefig('plots/GW150914_strain_whiten_filter.png')

#----------------------------------------------------------------
# Frequency analytic
#----------------------------------------------------------------

def gwfreq(iM,iT,iT0):
    const = (134.*np.pi*2)*np.power((1.21/iM),5./8.)
    output = const*np.power(np.maximum((iT0-iT),3e-2),-3./8.) # we can max it out above 500 Hz-ish
    return output

times = np.linspace(0., 4., 50)
freq = gwfreq(20, times, 4)
plt.figure()
plt.plot(times, freq)
plt.xlabel('time (s)')
plt.ylabel('Frequency (Hz)')
plt.title('Frequency variation in time of a black hole merger 20 $M_{\odot}$')
plt.savefig('plots/GW150914_spectrogram_analytic.png')

#----------------------------------------------------------------
# Spectrogram
#----------------------------------------------------------------

tevent = 1126259462.422         # Mon Sep 14 09:50:45 GMT 2015
indxt = np.where((sampletimes >= tevent-10) & (sampletimes < tevent+10))
NFFT = 4096/16
NOVL = NFFT*15/16
window = np.blackman(NFFT)
spec_cmap='viridis'

plt.figure()
spec_H1, freqs, bins, im = plt.specgram(white_data[indxt], NFFT=NFFT, Fs=4096, window=window, 
                                        noverlap=NOVL, cmap=spec_cmap, xextent=[-10,10])
plt.xlabel('time (s) since '+str(tevent))
plt.ylabel('Frequency (Hz)')
plt.colorbar()
plt.axis([-0.5, 0.5, 0, 500])
plt.title('aLIGO H1 strain data near GW150914')
plt.savefig('plots/GW150914_spectrogram_whitened.png')


#----------------------------------------------------------------
# Wave form analytic
#----------------------------------------------------------------

def osc(x,M,t0,n,phi):
    freq = gwfreq(M,x,t0)
    val = n*(np.cos(freq*(t0-x)+phi))*1e-12
    val = val*np.power(M*freq,10./3.)*(1*(x<=t0)+np.exp((freq/(2*np.pi))*(t0-x))*(x>t0))
    return val

def osc_dif(params, x, data, eps):
    iM=params["M"]
    iT0=params["t0"]
    norm=params["n"]
    phi=params["phi"]
    val=osc(x, iM, iT0, norm, phi)
    return (val-data)/eps

times = np.linspace(-0.1, 0.3, 1000)
freq = osc(times, 30, 0.18, 1, 0.0)
plt.figure()
plt.plot(times, freq)
plt.xlabel('Time (s)')
plt.ylabel('strain')
plt.title('Approximate GW waveform')
plt.savefig('plots/GW150914_strain_analytic.png')

#----------------------------------------------------------------
# Zoom and Fit
#----------------------------------------------------------------

indxt = np.where((sampletimes >= tevent-0.17) & (sampletimes < tevent+0.13))
x = sampletimes[indxt]
x = x-x[0]
white_data_zoom = white_data[indxt]
white_data_zoom_bp = white_data_bp[indxt]

plt.figure()
plt.plot(x, white_data_zoom_bp)
plt.xlabel('Time (s)')
plt.ylabel('strain')
plt.title('Advanced LIGO strain data near GW150914 filtered and zoomed')
plt.savefig('plots/GW150914_strain_whiten_filter_zoom.png')

import lmfit
from lmfit import Model,minimize, fit_report, Parameters

#Make a fit model using my favorite python fit lmfit
model = lmfit.Model(osc)
p = model.make_params()
p['M'].set(20)     # Mass guess
p['t0'].set(0.18)  # By construction we put the merger in the center
p['n'].set(1)      # normalization guess
p['phi'].set(0)    # Phase guess
unc = np.full(len(white_data_zoom),20)
out = minimize(osc_dif, params=p, args=(x, white_data_zoom, unc))
print(fit_report(out))

plt.plot(x, model.eval(params=out.params,x=x),'r',label='best fit')
plt.savefig('plots/GW150914_strain_whiten_filter_fit.png')

#----------------------------------------------------------------
# Significance vs. time
#----------------------------------------------------------------

def fitrange(data,xx,tcenter,trange):
    print(xx)
    print(tcenter, trange)
    findxt = np.where((xx >= tcenter-trange*0.5) & (xx < tcenter+trange*0.5))
    fwhite_data = data[findxt]
    x = xx[findxt]
    x = x-x[0]
    model = lmfit.Model(osc)
    p = model.make_params()
    p['M'].set(30)
    p['t0'].set(trange*0.5)
    p['n'].set(1)
    p['phi'].set(0)
    unc=np.full(len(fwhite_data),20)
    out = minimize(osc_dif, params=p, args=(x, fwhite_data, unc))
    return abs(out.params["n"].value/out.params["n"].stderr),out.redchi

times = np.arange(-14, 14, 0.05)
times += tevent
sigs=[]
chi2=[]
for time in times:
        pSig,pChi2 = fitrange(white_data_bp, sampletimes, time, 0.4)
        sigs.append(pSig)
        chi2.append(pChi2)

plt.figure()
plt.plot(times, sigs)
plt.xlabel('Time (s)')
plt.ylabel('N/$\sigma_{N}$')
plt.title('Significance of a potential wave')
plt.savefig('plots/GW150914_strain_significance.png')

plt.figure()
plt.plot(times, chi2)
plt.xlabel('Time (s)')
plt.ylabel('$\chi^{2}$')
plt.title('$\chi^{2}$ of a potential wave')
plt.savefig('plots/GW150914_strain_chi2.png')


