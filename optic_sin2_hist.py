import numpy as np
import pylab as pl
from pydelay import dde23

# define the equations
eqns = { 
        'x' : '(1/tau) *(-x-y+beta*pow(sin(x(t-T)+Phi0),2))', 
        'y' : '(1/theta)*x'
        }

bsteps=600;
Hist=np.zeros((bsteps,300));

for k in range(1,bsteps+1):

    #define the parameters, times is in 'ms'
    params = {
        'tau'   : 0.008,     # 0.008 ms < tau< 0.10 ms
        'beta'  : 0.5+k*0.005,
        'T'     : 10.5,       # 1.6 ms < T< 130 ms
        'Phi0'  : 0.25*np.pi,
        'theta' : 10.0         # 0.8 ms < theta< 1*1000 ms
        }
    #print(params)
    
    
    
    # Initialise the solver
    dde = dde23(eqns=eqns, params=params)
    
    #set the simulation parameters
    # (solve from t=0 to t=1000 and limit the maximum step size to 1.0)
    tfinal=6000
    tcut=700
    dde.set_sim_params(tfinal=tfinal, dtmax=0.5, AbsTol=10**-6, RelTol=10**-3)
    
    # set the history of to the constant function 0.5 (using a python lambda function)
    histfunc = {
        'x': lambda t: 0.3453,
        'y': lambda t: 0.82321
        }
    dde.hist_from_funcs(histfunc, 100)
    
    
    # run the simulator
    dde.run()
    
    T=params['T']
    beta=params['beta']
    sol1 = dde.sample((tfinal-tcut)+T, tfinal,0.1)
    x1 = sol1['x']
    H, xedges = np.histogram(x1, bins=300, range=(-2.0,2.0))
    Hist[k-1,:]=H.astype(float)/max(H.astype(float))

pl.imshow(np.transpose(Hist)) 
pl.axis('tight')
pl.show()

