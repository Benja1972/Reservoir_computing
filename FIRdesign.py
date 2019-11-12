import scipy as sp 
from pylab import *
import scipy.signal as signal
import numpy as np

#Plot frequency and phase response
def mfreqz(b,a=1):
    w,h = signal.freqz(b,a)
    h_dB = 20 * log10 (abs(h))
    subplot(211)
    plot(w/max(w),h_dB)
    ylim(-150, 5)
    ylabel('Magnitude (db)')
    xlabel(r'Normalized Frequency (x$\pi$ rad/sample)')
    title(r'Frequency response')
    subplot(212)
    h_Phase = unwrap(arctan2(imag(h),real(h)))
    plot(w/max(w),h_Phase)
    ylabel('Phase (radians)')
    xlabel(r'Normalized Frequency (x$\pi$ rad/sample)')
    title(r'Phase response')
    subplots_adjust(hspace=0.5)

#Plot step and impulse response
def impz(b,a=1):
    l = len(b)
    impulse = repeat(0.,l); impulse[0] =1.
    x = arange(0,l)
    response = signal.lfilter(b,a,impulse)
    subplot(211)
    stem(x, response)
    ylabel('Amplitude')
    xlabel(r'n (samples)')
    title(r'Impulse response')
    subplot(212)
    step = cumsum(response)
    stem(x, step)
    ylabel('Amplitude')
    xlabel(r'n (samples)')
    title(r'Step response')
    subplots_adjust(hspace=0.5)


# We want to construct a filter with a passband at 0.2-0.4 Hz, and
# stop bands at 0-0.1 Hz and 0.45-0.5 Hz. Note that this means that the
# behavior in the frequency ranges between those bands is unspecified and
# may overshoot.

bpass = signal.remez(72, [0, 0.1, 0.2, 0.4, 0.45, 0.5], [0, 1, 0],None,1./1)
freq, response = sp.signal.freqz(bpass)
ampl = np.abs(response)


fig = figure(100)
ax1 = fig.add_subplot(111)
ax1.semilogy(freq/(2*np.pi), ampl, 'b-') # freq in Hz
# [<matplotlib.lines.Line2D object at 0xf486790>]
#show()

#~ Lowpass FIR filter
#~ Designing a lowpass FIR filter is very simple to do with SciPy, 
#~ all you need to do is to define the window length, cut off frequency 
#~ and the window:

n = 100
alow = signal.firwin(n, cutoff = 0.8, window = "hamming")
#Frequency and phase response
figure(101)
mfreqz(alow)
#show()
#Impulse and step response
figure(102)
impz(alow)
#show()

#~ Highpass FIR Filter
#~ SciPy does not have a function for directly designing a highpass 
#~ FIR filter, however it is fairly easy design a lowpass filter and 
#~ use spectral inversion to convert it to highpass. See e.g Chp 16 of 
#~ The Scientist and Engineer's Guide to Digital Signal Processing for 
#~ the theory, the last page has an example code.


n = 101
a = signal.firwin(n, cutoff = 0.3, window = "hanning")
#Spectral inversion
a = -a
a[n/2] = a[n/2] + 1
figure(103)
mfreqz(a)
#show()


#~ Bandpass FIR filter
#~ To get a bandpass FIR filter with SciPy we first need to design 
#~ appropriate lowpass and highpass filters and then combine them:

n = 1001
#Lowpass filter
a = signal.firwin(n, cutoff = 0.3, window = 'blackmanharris')
#Highpass filter with spectral inversion
b = - signal.firwin(n, cutoff = 0.5, window = 'blackmanharris'); b[n/2] = b[n/2] + 1
#Combine into a bandpass filter
d = - (a+b); d[n/2] = d[n/2] + 1
#Frequency response
figure(104)
mfreqz(d)
show()


