from pylab import *
import scipy as sp 
import scipy.signal as signal
from scipy.stats import norm
import numpy as np


#Plot frequency and phase response

def nyqfreqz(b,nyq,a=1):
    w,h = signal.freqz(b,a)
    h_dB = 20 * log10 (abs(h))
    subplot(211)
    plot(nyq*w/max(w),h_dB)
    ylim(-150, 5)
    ylabel('Magnitude (db)')
    xlabel(r'Frequency (Hz)')
    title(r'Frequency response')
    subplot(212)
    h_Phase = unwrap(arctan2(imag(h),real(h)))
    plot(nyq*w/max(w),h_Phase)
    ylabel('Phase (radians)')
    xlabel(r'Frequency (Hz)')
    title(r'Phase response')
    subplots_adjust(hspace=0.5)

def FIRresponse(w,om,N=150):
    SE=np.zeros(np.size(om))
    for i in range(1,N+1):
        SE=SE+w[i-1]*np.exp(-1j*om*i*narmult*dtau)

    SE_dB = 20 * np.log10 (abs(SE))
    return SE_dB
    

N=150
theta = 1.59e-3
tau = 7.958e-6
left_band = 1/(2.*np.pi*theta)
right_band = 1/(2.*np.pi*tau)
dtau=0.2*tau
input_freq = 1/(2*np.pi*dtau)
nyq=1/(2.*dtau)
cutoffHz=20e3
narmult=1

f=np.arange(0.,0.1*nyq,1e-1)
om=2.*np.pi*f


#w = signal.firwin(N, cutoff = .95, window = "hamming")

#w1=np.random.rand(N,)
#w1=w1/np.sum(w1)
#w1=np.round(w1*2048.)


#w20 = 1.5*norm.rvs(size=N)+5.5
#w2=w20/np.sum(np.abs(w20))
#w2=np.round(w2*2048.)


w=sp.genfromtxt('w.dat')
#w=np.round(w*2048.)
#w1=ones(N,)



SE_dB = FIRresponse(w,om)
#SE_dB1 = FIRresponse(w1,om)
#SE_dB2 = FIRresponse(w2,om)

# Figures
plot_params = {'axes.labelsize': 18,
               'text.fontsize': 20,
               'legend.fontsize': 20,
               'title.fontsize': 22,
               'xtick.labelsize': 18,
               'ytick.labelsize': 18}
          
rcParams.update(plot_params)


#figure(90)
#plot(w20)


figure(100)
semilogx(f,SE_dB, 'b-')
#semilogx(f,SE_dB1, 'm-')
#semilogx(f,SE_dB2, 'g-')
semilogx([input_freq, input_freq],[-60,10],'m-',lw=2)
semilogx([left_band, left_band],[-60,10],'r-',lw=2)
semilogx([right_band, right_band],[-60,10],'g-',lw=2)
#print input_freq
#plot(f,SE_dB, 'r-')
grid(True, which='minor')
grid(True)
#axis('tight')
ylabel('Magnitude (db)')
xlabel(r'Frequency (Hz)')
title(r'Frequency response')


#~ figure(101)
#~ subplot(211)
#~ plot(f,SE_dB, 'b-')
#~ grid(True, which='minor')
#~ ylabel('Magnitude (db)')
#~ xlabel(r'Frequency (Hz)')
#~ title(r'Frequency response')
#~ 
#~ subplot(212)
#~ SE_Phase = unwrap(arctan2(imag(SE),real(SE)))
#~ plot(f,SE_Phase)
#~ ylabel('Phase (radians)')
#~ xlabel(r'Frequency (Hz)')
#~ title(r'Phase response')
#~ subplots_adjust(hspace=0.5)





#~ figure(103)
#~ plot(w)
#~ ylabel('Value in integer numbers')
#~ xlabel(r'Coefficients $w_i$')

#np.savetxt('w150(integer).dat', w)

show()


