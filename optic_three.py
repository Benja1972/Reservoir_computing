import numpy as np
import pylab as pl
#from math import sin, cos
from pydelay import dde23

# define the equations
#eqns = { 'x' : '1/tau * x(t-tau) / (1.0 + pow(x(t-tau),10.0)) -0.1*x' }

eqns = { 
        'x' : '(1/tau) *(-x-y+beta*pow(cos(x(t-T)-x(t-T-dT)+Phi0),2))', 
        'y' : '(1/theta)*x'
        }

#define the parameters
params = {
    'tau'   : 0.0122,
    'beta'  : 1.34191,
    'T'     : 24.35,
    'dT'    : 0.400,
    'Phi0'  : 0.23*np.pi,
    'theta' : 5.3*1000
    }

# Initialise the solver
dde = dde23(eqns=eqns, params=params)

#set the simulation parameters
# (solve from t=0 to t=1000 and limit the maximum step size to 1.0)
tfinal=17000
tcut=180
dde.set_sim_params(tfinal=tfinal, dtmax=0.5, AbsTol=10**-6, RelTol=10**-3)

# set the history of to the constant function 0.5 (using a python lambda function)
histfunc = {
    'x': lambda t: 0.3453,
    'y': lambda t: 0.0
    }
dde.hist_from_funcs(histfunc, 100)


# run the simulator
dde.run()

# Make a plot of x(t) vs x(t-tau):
# Sample the solution twice with a stepsize of dt=0.1:

# once in the interval
T=params['T']
dT=params['dT']
sol1 = dde.sample((tfinal-tcut)+T+dT, tfinal)
x1 = sol1['x']
y1 = sol1['y']
t = sol1['t']

# and once between
sol2 = dde.sample((tfinal-tcut)+dT, tfinal-T)
x2 = sol2['x']

# and once between
sol3 = dde.sample((tfinal-tcut), tfinal-T-dT)
x3 = sol3['x']
x4=x2-x3

pl.figure(1)
pl.subplot(311)
pl.plot(t,x1)
pl.xlabel('$t$')
pl.ylabel('$x(t)$')

pl.subplot(312)
pl.plot(t,y1)
pl.xlabel('$t$')
pl.ylabel('$y(t)$')


pl.subplot(313)
pl.plot(x4, x1)
pl.xlabel('$x(t-T)-x(t-T-\delta T)$')
pl.ylabel('$x(t)$')
#pl.show()

pl.figure(2)

H, xedges, yedges = np.histogram2d(x4, x1, bins=100)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
pl.imshow(H, extent=extent)
pl.show()



