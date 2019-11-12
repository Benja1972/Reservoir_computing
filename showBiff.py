import numpy as np
import pylab as pl

bn=0.1;
bf=8.0;
bsteps=200;
hb=(bf-bn)/bsteps;

params = {
    'tau'   : 1.0,                # 0.008 ms < tau< 0.10 ms
    'beta'  : bn,
    'T'     : 0.2*150,            # 1.6 ms < T< 130 ms
    'Phi0'  : 0.05,                    # 0.25*np.pi,
    'theta' : 603.2,              # 603.2,          ms
    'sc'    : 1/2048.0        # 1/72.3924
    }

T=params['T']
tau=params['tau']
theta=params['theta']


hbin=300;
# Hist=np.zeros((bsteps+1,hbin));
Hdwn=-2.0;
Hup=2.0;

nphi=12;
for j in range(1,nphi+1):
	#Phi0=np.pi*0.01+(j-1)*np.pi*((1.0-0.01)/nphi)
	#Hist = np.loadtxt('Biff_phi0 %0.3f.dat' %Phi0)
	
	Phi0=np.pi*0.005+(j-1)*np.pi*((0.5-0.005)/(nphi-1))
	Hist = np.loadtxt('Biff_phi0Scan01_%0.3f.dat' %Phi0)
	
	
	# Show picture
	pl.figure(int((j-1)/6)+1)
	pl.subplots_adjust(hspace=0.3)
	pl.suptitle(r'$\theta=$ %1.2f, $\tau=$ %1.2f, T= %1.2f' 
		% (theta, tau, T) )
	
	pl.subplot(3,2,(j-1)%6+1)
	pl.imshow(Hist, extent=(bn,bf,Hdwn,Hup)) 
	pl.axis('tight')
	pl.title(r'$\phi_0$= %1.2f $\pi$' 
	% (Phi0/np.pi) )
	pl.xlabel(r'$\beta$')
	pl.ylabel('$x(t)$')
	#pl.show()


pl.show()
