import numpy as np
import pylab as pl
from math import sin, cos
from pydelay import dde23

# define the equations
#eqns = { 'x' : '1/tau * x(t-tau) / (1.0 + pow(x(t-tau),10.0)) -0.1*x' }

eqns = { 'x' : '(1/tau) *(beta*pow(cos(x(t-T)-x(t-T-dT)+Phi0),2) -x)' }

#define the parameters
params = {
    'tau'   : 0.0122,
    'beta'  : 1.2,
    'T'     : 24.35,
    'dT'    : 0.400,
    'Phi0'  : 0.23*np.pi
    }

# Initialise the solver
dde = dde23(eqns=eqns, params=params)

#set the simulation parameters
# (solve from t=0 to t=1000 and limit the maximum step size to 1.0)
tfinal=10000
dde.set_sim_params(tfinal=tfinal, dtmax=0.5, AbsTol=10**-6, RelTol=10**-3)

# set the history of to the constant function 0.5 (using a python lambda function)
histfunc = {
    'x': lambda t: 0.3453
    }
dde.hist_from_funcs(histfunc, 100)


# run the simulator
dde.run()

# Make a plot of x(t) vs x(t-tau):
# Sample the solution twice with a stepsize of dt=0.1:

# once in the interval
T=params['T']
dT=params['dT']
sol1 = dde.sample((tfinal-80)+T+dT, tfinal, 0.01)
x1 = sol1['x']
t = sol1['t']

# and once between
sol2 = dde.sample((tfinal-80)+dT, tfinal-T, 0.01)
x2 = sol2['x']

# and once between
sol3 = dde.sample((tfinal-80), tfinal-T-dT, 0.01)
x3 = sol3['x']


pl.subplot(211)
pl.plot(t,x1)
pl.xlabel('$x(t)$')
pl.ylabel('$t$')


pl.subplot(212)
pl.plot(x2-x3, x1)
pl.xlabel('$x(t-T)-x(t-T-\delta T)$')
pl.ylabel('$x(t)$')

pl.show()
