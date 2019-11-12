import mdp

def extend(x,nl):
	y=mdp.numx.kron(mdp.numx.ones((nl,1)),x)
	y=y.T.reshape(1,-1)
	return y.T
	


class OpticDelayNode(mdp.Node):
		
	def __init__(self, input_dim=1, output_dim=None, dtype='float64', 
				  beta=0.8, gamma=0.4, phi0=3.14159265/4., theta=603.2, 
				  mask=None, w=None, neuron_width=0.2, bin=4):
		super(OpticDelayNode, self).__init__(input_dim=input_dim, output_dim=output_dim, dtype=dtype)
		
		self.beta = beta
		self.gamma = gamma
		self.phi0 = phi0
		self.theta = theta
		
		self.neuron_width = neuron_width
		self.bin = bin
		self.dt = float(self.neuron_width / self.bin)
		self.ntau = bin*self.output_dim
		self.ntrans = 500*self.ntau
		
		self.initial_state = mdp.numx.zeros((1, self.output_dim))
		self.x_hist_initial = 0.5*mdp.numx.random.random(( self.ntau, 1))
		self.y_hist_initial = 0.5*mdp.numx.random.random(( self.ntau, 1))
		
		self.mask_initial = mask
		self.w_initial = w
		
		self.mask = []
		self.w = []
		
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
			self.mask = (mdp.numx.ceil(mdp.numx.random.rand(self.output_dim, self.input_dim)/0.5) * 2 - 3)
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

	def _get_supported_dtypes(self):
		return ['float32', 'float64']
    
	def _execute(self, x):
		steps = x.shape[0]
		
		# Pre-allocate the state vector, adding the initial state
		states = mdp.numx.concatenate((self.initial_state, mdp.numx.zeros((steps, self.output_dim))))
		
		self.x_hist = self.x_hist_initial.copy()
		self.y_hist = self.y_hist_initial.copy()
		
		# Loop over the input data and compute the reservoir states
		for n in range(steps):
			inp=extend(self.gamma*mdp.numx.dot(self.mask,x[n,:]),self.bin)
			xt=mdp.numx.zeros(( self.ntau+self.ntau,1))
			yt=mdp.numx.zeros(( self.ntau+self.ntau,1))
		
			xt[:self.ntau,:]=self.x_hist.copy()
			yt[:self.ntau,:]=self.y_hist.copy()
			
			for i in range(self.ntau-1,self.ntau+self.ntau-1):
				xd=mdp.numx.dot(self.w,xt[i-self.ntau+1:i+1:self.bin,:])
				
				xt[i+1,:]=xt[i,:]+self.dt*(-xt[i,:]-(1/self.theta)*yt[i,:]+\
				self.beta*mdp.numx.sin(xd+inp[i-self.ntau+1,:]+self.phi0)**2)
				
				yt[i+1,:]=yt[i,:]+self.dt*xt[i,:]
			
			states[n + 1, :] = xt[-self.ntau+int(self.bin/2)::self.bin,:].copy().T
			
			self.x_hist = xt[-self.ntau:,:].copy()
			self.y_hist = yt[-self.ntau:,:].copy()
			self.inp = inp

		return states[1:, :]




class TimeTwoNode(mdp.Node):
	
	def is_trainable(self):
		return False
		
	def _execute(self,x):
		return 2*x
		
	def _inverse(self,x):
		return y/2
		
	
class PowerNode(mdp.Node):
	def __init__(self, power, input_dim=None, dtype=None):
		super(PowerNode, self).__init__(input_dim=input_dim, dtype=dtype)
		self.power = power
		
	def is_trainable(self):
		return False
		
	def is_invertible(self):
		return False
		
	def _get_supported_dtypes(self):
		return ['float32', 'float64']
		
	def _execute(self, x):
		return self._refcast(x**self.power)

class UnitVarianceNode(mdp.Node):
     def __init__(self, input_dim=None, dtype=None):
         super(UnitVarianceNode, self).__init__(input_dim=input_dim,
                                                 dtype=dtype)
         self.avg = None # average
         self.std = None # standard deviation
         self.tlen = 0
     def _get_train_seq(self):
         return [(self._train_mean, self._stop_mean),
                 (self._train_std, self._stop_std)]
     def _train_mean(self, x):
         if self.avg is None:
             self.avg = mdp.numx.zeros(self.input_dim,
                                       dtype=self.dtype)
         self.avg += mdp.numx.sum(x, 0)
         self.tlen += x.shape[0]
     def _stop_mean(self):
         self.avg /= self.tlen
     def _train_std(self, x):
         if self.std is None:
             self.tlen = 0
             self.std = mdp.numx.zeros(self.input_dim,
                                       dtype=self.dtype)
         self.std += mdp.numx.sum((x - self.avg)**2., 0)
         self.tlen += x.shape[0]
     def _stop_std(self):
         # compute the standard deviation
         self.std = mdp.numx.sqrt(self.std/(self.tlen-1))
     def _execute(self, x):
         return (x - self.avg)/self.std
     def _inverse(self, y):
         return y*self.std + self.avg
