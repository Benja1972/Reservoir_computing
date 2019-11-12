import numpy as np
import pylab as pl
from pydelay import dde23

# define the equations
eqns = { 
        'x' : '(1/tau) *(-x-y+beta*pow(sin(x(t-T)+Phi0),2))',
        'y' : '(1/theta)*x'
        }

# total delay
T = 10.5
dt = T/1000
nsteps =10 
eps_n = 0.2*np.random.random() 

params = {
    'tau'   : 0.08,     # 0.008 ms < tau< 0.10 ms
    'beta'  : 2.4,
    'T'     : T,       # 1.6 ms < T< 130 ms
    'Phi0'  : 0.25*np.pi,
    'theta' : 543.0,         # 0.8 ms < theta< 1*1000 ms
    }
print(params)



thist = np.linspace(-T, 0, 200)
xhist = 0.3453*np.ones(len(thist))
yhist = 0.82321*np.ones(len(thist))

hist = {
    't' : thist,
    'x' : xhist,
    'y' : yhist
}


# Initialise the solver
dde = dde23(eqns=eqns, params=params)

tfinal=nsteps*T
tcut=(nsteps)*T
#set the simulation parameters
dde.set_sim_params(tfinal=tfinal, dtmax=0.15, AbsTol=10**-6, RelTol=10**-3)

# set the history 
dde.hist_from_arrays(hist)

# run the simulator
dde.run()

sol = dde.sample((tfinal-tcut)+T, tfinal,0.1)
x = sol['x']
y = sol['y']
t = sol['t']


    


# Make a plot of x(t) vs x(t-tau):
# Sample the solution twice with a stepsize of dt=0.1:


#~ T=params['T']
#~ beta=params['beta']
#~ tau=params['tau']
#~ theta=params['theta']


#~ sol1 = dde.sample((tfinal-tcut)+T, tfinal,0.1)
#~ x1 = sol1['x']
#~ y1 = sol1['y']
#~ t = sol1['t']
#~ 
#~ # and once between
#~ sol2 = dde.sample((tfinal-tcut), tfinal-T,0.1)
#~ x2 = sol2['x']


pl.figure(1)


pl.subplot(211)
pl.plot(t,x)
pl.xlabel('$t$')
pl.ylabel('$x(t)$')
pl.title(r'$\beta=$ %1.2f,  $\theta=$ %1.2f, $\tau=$ %1.2f, T= %1.2f' 
            % (params['beta'], params['theta'], params['tau'], T ) )

pl.subplot(212)
pl.plot(t,y)
pl.xlabel('$t$')
pl.ylabel('$y(t)$')

pl.show()




