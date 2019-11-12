import numpy as np
import pylab as pl
from math import sin, cos
from pydelay import dde23

# define the equations
eqns = { 
        'x' : '-x-y+gamma*(pow(cos(m+d*tanh(x(t-tau))),2)-pow(cos(m),2))',
        'y' : 'eps*x'
        }

#noise = { 'x': 'ampn * gwn()'}

#define the parameters, times is in 'ms'
params = {
    'tau'   : 1820.0,       
    'gamma' : 4.52,       # 0.00 < gamma< 5.00
    'm'     : 0.0*np.pi,  # -pi/2 < m < pi/2 
    'd'     : 2.1,     
    'eps'   : 2.0*10**-6,  
    }
print(params)

# Initialise the solver
dde = dde23(eqns=eqns, params=params)

#set the simulation parameters
# (solve from t=0 to t=1000 and limit the maximum step size to 1.0)
tfinal=70000
tcut=1000
dde.set_sim_params(tfinal=tfinal, dtmax=1.0, AbsTol=10**-6, RelTol=10**-3)

# set the history of to the constant function 0.5 (using a python lambda function)
histfunc = {
    'x': lambda t: 0.2*sin(.2*t),
    'y': lambda t: 0.2*cos(0.00001*t)
    }
dde.hist_from_funcs(histfunc, 40000)


# run the simulator
dde.run()

# Make a plot of x(t) vs x(t-tau):
# Sample the solution twice with a stepsize of dt=0.1:


tau=params['tau']
gamma=params['gamma']
m=params['m']

x=dde.sol['x']
t=dde.sol['t']


pl.figure(1)

pl.plot(t,x)

pl.show()



#sol1 = dde.sample((tfinal-tcut)+tau, tfinal,0.1)
#x1 = sol1['x']
#y1 = sol1['y']
#t = sol1['t']
#
## and once between
#sol2 = dde.sample((tfinal-tcut), tfinal-tau,0.1)
#x2 = sol2['x']
#
#
#pl.figure(1)
#
#pl.subplot(311)
#pl.plot(t,x1)
#pl.xlabel('$t$')
#pl.ylabel('$x(t)$')
#pl.title(r'$\gamma=$ %1.2f,  $\tau=$ %1.2f, $m=$ %1.2f' 
#            % (gamma, tau, m ) )
#
#pl.subplot(312)
#pl.plot(t,y1)
#pl.xlabel('$t$')
#pl.ylabel('$y(t)$')
#
#
#pl.subplot(313)
#pl.plot(x2, x1,'.')
#pl.xlabel(r'$x(t-\tau)$')
#pl.ylabel('$x(t)$')
#
#pl.figure(2)
#
#H, xedges, yedges = np.histogram2d(x1, x2, bins=100)
#extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#pl.imshow(H, extent=extent)
#pl.xlabel('$x(t)$')
#pl.ylabel(r'$x(t-\tau)$')
#
#pl.show()
#
#

