import numpy as np
import pylab as pl
from pydelay import dde23

# define the equations
eqns = { 
        'x' : '(-x-(1/theta)*y+beta*pow(sin(sc*(16.00000*x(t-0.0067*T)+\
                                                6.00000*x(t-0.0133*T)+\
                                                22.00000*x(t-0.0200*T)+\
                                                25.00000*x(t-0.0267*T)+\
                                                22.00000*x(t-0.0333*T)+\
                                                13.00000*x(t-0.0400*T)+\
                                                7.00000*x(t-0.0467*T)+\
                                                19.00000*x(t-0.0533*T)+\
                                                6.00000*x(t-0.0600*T)+\
                                                25.00000*x(t-0.0667*T)+\
                                                16.00000*x(t-0.0733*T)+\
                                                15.00000*x(t-0.0800*T)+\
                                                4.00000*x(t-0.0867*T)+\
                                                2.00000*x(t-0.0933*T)+\
                                                7.00000*x(t-0.1000*T)+\
                                                22.00000*x(t-0.1067*T)+\
                                                23.00000*x(t-0.1133*T)+\
                                                18.00000*x(t-0.1200*T)+\
                                                19.00000*x(t-0.1267*T)+\
                                                6.00000*x(t-0.1333*T)+\
                                                15.00000*x(t-0.1400*T)+\
                                                21.00000*x(t-0.1467*T)+\
                                                10.00000*x(t-0.1533*T)+\
                                                25.00000*x(t-0.1600*T)+\
                                                2.00000*x(t-0.1667*T)+\
                                                8.00000*x(t-0.1733*T)+\
                                                13.00000*x(t-0.1800*T)+\
                                                2.00000*x(t-0.1867*T)+\
                                                19.00000*x(t-0.1933*T)+\
                                                14.00000*x(t-0.2000*T)+\
                                                14.00000*x(t-0.2067*T)+\
                                                21.00000*x(t-0.2133*T)+\
                                                22.00000*x(t-0.2200*T)+\
                                                20.00000*x(t-0.2267*T)+\
                                                8.00000*x(t-0.2333*T)+\
                                                12.00000*x(t-0.2400*T)+\
                                                19.00000*x(t-0.2467*T)+\
                                                3.00000*x(t-0.2533*T)+\
                                                3.00000*x(t-0.2600*T)+\
                                                7.00000*x(t-0.2667*T)+\
                                                13.00000*x(t-0.2733*T)+\
                                                25.00000*x(t-0.2800*T)+\
                                                18.00000*x(t-0.2867*T)+\
                                                8.00000*x(t-0.2933*T)+\
                                                7.00000*x(t-0.3000*T)+\
                                                22.00000*x(t-0.3067*T)+\
                                                23.00000*x(t-0.3133*T)+\
                                                16.00000*x(t-0.3200*T)+\
                                                7.00000*x(t-0.3267*T)+\
                                                2.00000*x(t-0.3333*T)+\
                                                22.00000*x(t-0.3400*T)+\
                                                15.00000*x(t-0.3467*T)+\
                                                24.00000*x(t-0.3533*T)+\
                                                2.00000*x(t-0.3600*T)+\
                                                15.00000*x(t-0.3667*T)+\
                                                7.00000*x(t-0.3733*T)+\
                                                21.00000*x(t-0.3800*T)+\
                                                5.00000*x(t-0.3867*T)+\
                                                11.00000*x(t-0.3933*T)+\
                                                10.00000*x(t-0.4000*T)+\
                                                21.00000*x(t-0.4067*T)+\
                                                17.00000*x(t-0.4133*T)+\
                                                5.00000*x(t-0.4200*T)+\
                                                8.00000*x(t-0.4267*T)+\
                                                3.00000*x(t-0.4333*T)+\
                                                17.00000*x(t-0.4400*T)+\
                                                15.00000*x(t-0.4467*T)+\
                                                4.00000*x(t-0.4533*T)+\
                                                4.00000*x(t-0.4600*T)+\
                                                12.00000*x(t-0.4667*T)+\
                                                23.00000*x(t-0.4733*T)+\
                                                14.00000*x(t-0.4800*T)+\
                                                1.00000*x(t-0.4867*T)+\
                                                1.00000*x(t-0.4933*T)+\
                                                21.00000*x(t-0.5000*T)+\
                                                12.00000*x(t-0.5067*T)+\
                                                10.00000*x(t-0.5133*T)+\
                                                20.00000*x(t-0.5200*T)+\
                                                9.00000*x(t-0.5267*T)+\
                                                14.00000*x(t-0.5333*T)+\
                                                18.00000*x(t-0.5400*T)+\
                                                22.00000*x(t-0.5467*T)+\
                                                8.00000*x(t-0.5533*T)+\
                                                17.00000*x(t-0.5600*T)+\
                                                25.00000*x(t-0.5667*T)+\
                                                2.00000*x(t-0.5733*T)+\
                                                15.00000*x(t-0.5800*T)+\
                                                11.00000*x(t-0.5867*T)+\
                                                8.00000*x(t-0.5933*T)+\
                                                7.00000*x(t-0.6000*T)+\
                                                19.00000*x(t-0.6067*T)+\
                                                26.00000*x(t-0.6133*T)+\
                                                5.00000*x(t-0.6200*T)+\
                                                20.00000*x(t-0.6267*T)+\
                                                5.00000*x(t-0.6333*T)+\
                                                25.00000*x(t-0.6400*T)+\
                                                21.00000*x(t-0.6467*T)+\
                                                11.00000*x(t-0.6533*T)+\
                                                19.00000*x(t-0.6600*T)+\
                                                13.00000*x(t-0.6667*T)+\
                                                21.00000*x(t-0.6733*T)+\
                                                9.00000*x(t-0.6800*T)+\
                                                2.00000*x(t-0.6867*T)+\
                                                15.00000*x(t-0.6933*T)+\
                                                23.00000*x(t-0.7000*T)+\
                                                5.00000*x(t-0.7067*T)+\
                                                11.00000*x(t-0.7133*T)+\
                                                19.00000*x(t-0.7200*T)+\
                                                1.00000*x(t-0.7267*T)+\
                                                24.00000*x(t-0.7333*T)+\
                                                20.00000*x(t-0.7400*T)+\
                                                14.00000*x(t-0.7467*T)+\
                                                5.00000*x(t-0.7533*T)+\
                                                13.00000*x(t-0.7600*T)+\
                                                13.00000*x(t-0.7667*T)+\
                                                26.00000*x(t-0.7733*T)+\
                                                22.00000*x(t-0.7800*T)+\
                                                25.00000*x(t-0.7867*T)+\
                                                17.00000*x(t-0.7933*T)+\
                                                10.00000*x(t-0.8000*T)+\
                                                24.00000*x(t-0.8067*T)+\
                                                12.00000*x(t-0.8133*T)+\
                                                6.00000*x(t-0.8200*T)+\
                                                10.00000*x(t-0.8267*T)+\
                                                18.00000*x(t-0.8333*T)+\
                                                14.00000*x(t-0.8400*T)+\
                                                19.00000*x(t-0.8467*T)+\
                                                26.00000*x(t-0.8533*T)+\
                                                25.00000*x(t-0.8600*T)+\
                                                14.00000*x(t-0.8667*T)+\
                                                25.00000*x(t-0.8733*T)+\
                                                3.00000*x(t-0.8800*T)+\
                                                1.00000*x(t-0.8867*T)+\
                                                8.00000*x(t-0.8933*T)+\
                                                15.00000*x(t-0.9000*T)+\
                                                14.00000*x(t-0.9067*T)+\
                                                23.00000*x(t-0.9133*T)+\
                                                14.00000*x(t-0.9200*T)+\
                                                11.00000*x(t-0.9267*T)+\
                                                14.00000*x(t-0.9333*T)+\
                                                18.00000*x(t-0.9400*T)+\
                                                0.00000*x(t-0.9467*T)+\
                                                21.00000*x(t-0.9533*T)+\
                                                4.00000*x(t-0.9600*T)+\
                                                12.00000*x(t-0.9667*T)+\
                                                7.00000*x(t-0.9733*T)+\
                                                9.00000*x(t-0.9800*T)+\
                                                17.00000*x(t-0.9867*T)+\
                                                4.00000*x(t-0.9933*T)+\
                                                7.00000*x(t-1.0000*T))+Phi0),2))',
        'y' : 'x'
        }

bn=0.1;
bf=8.0;
bsteps=200;
hb=(bf-bn)/bsteps;

# define the parameters, times is in 'ms'
params = {
    'tau'   : 1.0,                # 0.008 ms < tau< 0.10 ms
    'beta'  : bn,
    'T'     : 0.2*150,            # 1.6 ms < T< 130 ms
    'Phi0'  : 0.05,                    # 0.25*np.pi,
    'theta' : 603.2,              # 603.2,          ms
    'sc'    : 1/2048.0        # 1/72.3924
    }
# print(params)

# set the history of to the constant function 0.5 (using a python lambda function)
histfunc = {
        'x': lambda t: 0.2*np.sin(0.2*t),
        'y': lambda t: 0.1*np.cos(1.2*t)
    }


#output diagramm
hbin=300;
Hist=np.zeros((bsteps+1,hbin));
Hdwn=-2.0;
Hup=2.0;

nphi=12;
for j in range(1,nphi+1):
	params['Phi0']=np.pi*0.005+(j-1)*np.pi*((0.5-0.005)/(nphi-1))

	for k in range(1,bsteps+2):

		params['beta']=bn+(k-1)*hb
			
		# Initialise the solver
		dde = dde23(eqns=eqns, params=params)
		
		#set the simulation parameters
		# (solve from t=0 to t=tfinal and limit the maximum step size to dtmax)
		tfinal=20000
		tcut=10000
		dde.set_sim_params(tfinal=tfinal, dtmax=0.5, AbsTol=10**-6, RelTol=10**-3)
		
		# set the history
		dde.hist_from_funcs(histfunc, 100)
		
		
		# run the simulator
		dde.run()
	   
		x = dde.sol['x']
		t = dde.sol['t']
		H, xedges = np.histogram(x[np.ceil(len(x)*0.6):], bins=hbin, range=(Hdwn,Hup))
		Hist[k-1,:]=H.astype(float)/max(H.astype(float))


	# Save to file
	Phi0=params['Phi0']
	np.savetxt('Biff_phi0Scan01_%0.3f.dat' %Phi0 , Hist.T)

# Show picture
T=params['T']
#beta=params['beta']
tau=params['tau']
theta=params['theta']

pl.imshow(Hist.T, extent=(bn,bf,Hdwn,Hup)) 
pl.axis('tight')
pl.title(r'$\theta=$ %1.2f, $\tau=$ %1.2f, T= %1.2f, $\phi_0$= %1.2f' 
            % (theta, tau, T, Phi0) )
pl.xlabel(r'$\beta$')
pl.ylabel('$x(t)$')
pl.show()