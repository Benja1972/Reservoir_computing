import mdp
from scipy import weave
from scipy.weave import converters
import numpy as np
#import psyco
#psyco.full()

def extend(x,nl):
    y=mdp.numx.kron(mdp.numx.ones((nl,1)),x)
    y=y.T.reshape(1,-1)
    return y.T
    


class OpticDelayNode(mdp.Node):
        
    def __init__(self, input_dim=1, output_dim=None, dtype='float64', _instance=0,
                  beta=0.5, gamma=0.4, phi0=mdp.numx.pi/4., theta=603.2, 
                  mask=None, w=None, neuron_width=0.2, bin=4):
        super(OpticDelayNode, self).__init__(input_dim=input_dim, output_dim=output_dim, dtype=dtype)
        
        self.beta = beta
        self.gamma = gamma
        self.phi0 = phi0
        self.theta = theta
        
        self.neuron_width = neuron_width
        self.bin = bin
        self.dt = self.neuron_width / self.bin
        self.ntau = bin*self.output_dim
        self.ntrans = 500*self.ntau
        
        self.initial_state = mdp.numx.zeros((1, self.output_dim))
        self.x_hist_initial = 0.000*mdp.numx.random.random(( self.ntau, 1))
        self.y_hist_initial = 0.000*mdp.numx.random.random(( self.ntau, 1))
        
        self.mask_initial = mask
        self.w_initial = w
        
        self.mask = []
        self.w = []
        
        # Instance ID, used for making different reservoir instantiations with the same parameters
        # Can be ranged over to simulate different 'runs'
        self._instance = _instance
        
        self._is_initialized = False

        if input_dim is not None:
        # Call the initialize function to create the weight matrices
            self.initialize()

    def is_trainable(self):
        return False

    def is_invertible(self):
        return False
    
    def initialize(self):
        # Initialize input weight matrix(mask)
        if self.mask_initial is None:
            # Initialize it to uniform random values using input_scaling
            # self.mask = (mdp.numx.random.rand(self.output_dim, self.input_dim) * 2 - 1)
            self.mask = 0.1*(mdp.numx.ceil(mdp.numx.random.rand(self.output_dim, self.input_dim)/0.5) * 2 - 3)
        else:
            self.mask = self.mask_initial # else just copy it
    
        # Initialize reservoir weight matrix
        if self.w_initial is None:
            self.w = (mdp.numx.random.rand(1, self.output_dim) * 2 - 1)
            # scale it to spectral radius
            self.w /= mdp.numx.sum(self.w, axis=1)
        else:
            self.w = self.w_initial   # else just copy it
        
        xt=mdp.numx.zeros(( self.ntrans+self.ntau,1))
        yt=mdp.numx.zeros(( self.ntrans+self.ntau,1))
        
        xt[:self.ntau,:]=self.x_hist_initial.copy()
        yt[:self.ntau,:]=self.y_hist_initial.copy()
        
        for i in range(self.ntau-1,self.ntau+self.ntrans-1):
            xd=mdp.numx.dot(self.w,xt[i-self.ntau+1:i+1:self.bin,:])
            xt[i+1,:]=xt[i,:]+self.dt*(-xt[i,:]-(1/self.theta)*yt[i,:]+self.beta*mdp.numx.sin(xd+self.phi0)**2)
            yt[i+1,:]=yt[i,:]+self.dt*xt[i,:]
            
        self.x_hist_initial = xt[-self.ntau:,:].copy()
        self.y_hist_initial = yt[-self.ntau:,:].copy()
        self.xt = xt;
        self.yt = yt;
        
        self._is_initialized = True

    def _get_supported_dtypes(self):
        return ['float32', 'float64']
    
    def _execute(self, x):
        """ Executes simulation with input vector x.
        """
        # Check if the weight matrices are intialized, otherwise create them
        if not self._is_initialized:
            self.initialize()
        steps = x.shape[0]
        
        # Pre-allocate the state vector, adding the initial state
        states = mdp.numx.concatenate((self.initial_state, mdp.numx.zeros((steps, self.output_dim))))
        
        self.x_hist = self.x_hist_initial.copy()
        self.y_hist = self.y_hist_initial.copy()
        
        # Loop over the input data and compute the reservoir states
        for n in range(steps):
            inp=self.gamma*extend(mdp.numx.dot(self.mask,x[n,:]),self.bin)
            xt=mdp.numx.zeros(( self.ntau+self.ntau,1))
            yt=mdp.numx.zeros(( self.ntau+self.ntau,1))
        
            xt[:self.ntau,:]=self.x_hist.copy()
            yt[:self.ntau,:]=self.y_hist.copy()
            
            for i in range(self.ntau-1,self.ntau+self.ntau-1):
                xd=mdp.numx.dot(self.w,xt[i-self.ntau+1:i+1:self.bin,:])
                
                xt[i+1,:]=xt[i,:]+self.dt*(-xt[i,:]-(1/self.theta)*yt[i,:]+\
                self.beta*mdp.numx.sin(xd+inp[i-self.ntau+1,:]+self.phi0)**2)
                
                yt[i+1,:]=yt[i,:]+self.dt*xt[i,:]
            
            states[n + 1, :] = xt[-self.ntau+int(self.bin/2)+1::self.bin,:].copy().T
            #states[n , :] = xt[-self.ntau+int(self.bin/2)+1::self.bin,:].copy().T
            
            self.x_hist = xt[-self.ntau:,:].copy()
            self.y_hist = yt[-self.ntau:,:].copy()
            self.inp = inp

        return states[1:, :]
