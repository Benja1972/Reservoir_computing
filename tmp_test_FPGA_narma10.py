import Oger
import mdp
import pylab
import delay_nodes
import numpy as np

if __name__ == "__main__":

    inputs = 1
    neurons = 150
    n_samples = 4
    
    # parameters of model
    beta = 0.683
    gamma = 0.073
    phi0 = np.pi/4. 
    w = np.array(np.loadtxt('w.dat'), ndmin=2)
    
    [inputs, outputs] = Oger.datasets.narma10(n_samples=n_samples, sample_len=2*neurons)
    data = [inputs, zip(inputs, outputs)]
    #data = [ inputs, outputs]

    
    # construct individual nodes
    reservoir = delay_nodes.OpticDelayNode(output_dim=neurons, beta=beta, gamma=gamma, phi0=phi0, w=w)
    #readout = Oger.nodes.RidgeRegressionNode(ridge_param=0.000000000000001)
    #readout = mdp.nodes.LinearRegressionNode(use_pinv=True)
    #readout = mdp.nodes.RidgeScikitsLearnNode(alpha=0.001) 

    # build network with MDP framework
    #flow = Oger.nodes.InspectableFlow([reservoir, readout], verbose=1)
    flow = reservoir + readout
      
    # train the flow 
    flow.train(data)

    #apply the trained flow to the training data and test data
    trainout = flow(inputs[0])
    testout = flow(inputs[3])

    print "NRMSE: " + str(Oger.utils.nrmse(outputs[3], testout))

    #plot the input
    pylab.subplot(3, 1, 1)
    pylab.plot(inputs[0])
    
    pylab.subplot(3, 1, 2)
    pylab.plot(trainout, 'r')
    pylab.plot(outputs[0], 'b')

    pylab.subplot(3, 1, 3)
    pylab.plot(testout, 'r')
    pylab.plot(outputs[3], 'b')
    pylab.show()
    
    
    #~ pylab.subplot(nx, ny, 4)
    #~ pylab.plot(flow.inspect(reservoir))
    #~ pylab.show()
    #~ 
