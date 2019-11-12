import numpy as np
import pylab as pl
from pydelay import dde23

# define the equations
eqns = { 
        'x' : '(-x-(1/theta)*y+beta*pow(sin(x(t-T)+gamma*inp(t,delt,T)+Phi0),2))',
        'y' : 'x'
        }

#define the parameters, times is in 'ms'
params = {
    'tau'   : 1.0,          # 0.008 ms < tau< 0.10 ms
    'beta'  : .687000,
    'T'     : 0.2*150,          # 1.6 ms < T< 130 ms
    'delt'  : 0.1,
    'gamma' : 1.0,
    'Phi0'  : 0.25*np.pi,
    'theta' : 603.2,  #603.2,         # 0.8 ms < theta< 1*1000 ms
    }

inputcode = """
double inp(double t, double delt, double T) {
int n, c;
double  narma;
char line[80];
bool work;
FILE *fr;

work = true;
n = (int) ((t-10*T)/delt);
c= 0;
 fr = fopen("input.dat", "rt");  

   while(fgets(line, 80, fr) != NULL && work)
   {
         sscanf (line, "%lf", &narma);
         if( c==n ){
               work = false;
           }else{ c++; }
 }
fclose(fr);
return narma;
}
"""


print(params)

# Initialise the solver
dde = dde23(eqns=eqns, params=params, supportcode=inputcode)

#set the simulation parameters
# (solve from t=0 to t=tfinal and limit the maximum step size to dtmax)
T=params['T']
tfinal=50*T
tcut=47*T
dde.set_sim_params(tfinal=tfinal, dtmax=0.5, AbsTol=10**-1, RelTol=10**-1)

# set the history  using a python lambda function
histfunc = {
    'x': lambda t: 0.00,
    'y': lambda t: 207.20
    }
dde.hist_from_funcs(histfunc, 1000)


# run the simulator
dde.run()


T=params['T']
beta=params['beta']
tau=params['tau']
theta=params['theta']

# Solution and delayed solution
sol1 = dde.sample((tfinal-tcut)+T, tfinal,0.1)
x1 = sol1['x']
y1 = sol1['y']
t = sol1['t']


sol2 = dde.sample((tfinal-tcut), tfinal-T,0.1)
x2 = sol2['x']


# Figures
plot_params = {'axes.labelsize': 18,
               'font.size': 20,
               'legend.fontsize': 20,
               #~ 'title.fontsize': 22,
               'xtick.labelsize': 18,
               'ytick.labelsize': 18}
          
pl.rcParams.update(plot_params)

pl.figure(1)

pl.subplot(211)
pl.plot(t,x1)
pl.xlabel('$t$')
pl.ylabel('$x(t)$')
pl.title(r'$\beta=$ %1.3f,  $\theta=$ %1.2f, $\tau=$ %1.2f, T= %1.2f' 
            % (beta, theta, tau, T ) )

pl.subplot(212)
pl.plot(t,y1)
pl.xlabel('$t$')
pl.ylabel('$y(t)$')


#pl.subplot(313)
#pl.plot(x2, x1,'.')
#pl.xlabel('$x(t-T)$')
#pl.ylabel('$x(t)$')

#pl.figure(2)
#
#H, xedges, yedges = np.histogram2d(x1, x2, bins=100)
#extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#pl.imshow(H, extent=extent)
#pl.xlabel('$x(t)$')
#pl.ylabel('$x(t-T)$')

pl.show()


