import numpy as np
import pylab as pl
from pydelay import dde23

# define the equations
eqns = { 
        'x' : '(1/tau) *(-x-y+beta*pow(sin(x(t-T)+in+Phi0),2))',
        'y' : '(1/theta)*x'
        }

# total delay
T = 0.2*150
dt=T/1000
nsteps = 200

thist = np.linspace(-T, 0, 200)
xhist = 0.3453*np.ones(len(thist))
yhist = 0.82321*np.ones(len(thist))

hist = {
    't' : thist,
    'x' : xhist,
    'y' : yhist
}

tsol = thist
xsol = xhist
ysol = yhist

for i in range(nsteps):
    gamma = 0.00000000*(i<100)+0.1*(i>=100)
    input = gamma*(np.random.random()-0.5)*2. 
    #define the parameters, times is in 'ms'
    params = {
        'tau'   : 1.0,     # 0.008 ms < tau< 0.10 ms
        'beta'  : 0.68,
        'T'     : T,       # 1.6 ms < T< 130 ms
        'Phi0'  : 0.25*np.pi,
        'theta' : 603.2,         # 0.8 ms < theta< 1*1000 ms
        'in'    : input
        }
    print(params)

    # Initialise the solver
    dde = dde23(eqns=eqns, params=params)

    #set the simulation parameters
    dde.set_sim_params(tfinal=T, dtmax=0.15, AbsTol=10**-6, RelTol=10**-3)

    # set the history of to the constant function 0.5 (using a python lambda function)
    dde.hist_from_arrays(hist)

    # run the simulator
    dde.run()
    sol = dde.sample(tfinal=T+dt, dt=dt)
    hist = sol
    tsol = np.hstack((tsol,sol['t']+i*T))
    xsol = np.hstack((xsol,sol['x']))
    ysol = np.hstack((ysol,sol['y']))
    
    


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
pl.plot(tsol,xsol)
pl.xlabel('$t$')
pl.ylabel('$x(t)$')
pl.title(r'$\beta=$ %1.2f,  $\theta=$ %1.2f, $\tau=$ %1.2f, T= %1.2f' 
            % (params['beta'], params['theta'], params['tau'], T ) )

pl.subplot(212)
pl.plot(tsol,ysol)
pl.xlabel('$t$')
pl.ylabel('$y(t)$')

pl.show()

#~ pl.subplot(311)
#~ pl.plot(t,x1)
#~ pl.xlabel('$t$')
#~ pl.ylabel('$x(t)$')
#~ pl.title(r'$\beta=$ %1.2f,  $\theta=$ %1.2f, $\tau=$ %1.2f, T= %1.2f' 
            #~ % (beta, theta, tau, T ) )
#~ 
#~ pl.subplot(312)
#~ pl.plot(t,y1)
#~ pl.xlabel('$t$')
#~ pl.ylabel('$y(t)$')
#~ 
#~ 
#~ pl.subplot(313)
#~ pl.plot(x2, x1,'.')
#~ pl.xlabel('$x(t-T)$')
#~ pl.ylabel('$x(t)$')
#~ 
#~ pl.figure(2)
#~ 
#~ H, xedges, yedges = np.histogram2d(x1, x2, bins=100)
#~ extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#~ pl.imshow(H, extent=extent)
#~ pl.xlabel('$x(t)$')
#~ pl.ylabel('$x(t-T)$')

#~ pl.show()



