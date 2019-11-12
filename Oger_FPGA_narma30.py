import Oger
import pylab
import delay_nodes
import numpy as np

if __name__ == "__main__":
    """ This example demonstrates a very simple reservoir+readout setup on the 30th order NARMA task.
    """
    inputs = 1
    neurons = 150
    w = np.array(np.loadtxt('w.dat'), ndmin=2)
    #timesteps = 10000
    #washout = 30

    nx = 4
    ny = 1

    #[x, y] = Oger.datasets.narma30(sample_len=1000)
    [x, y] = Oger.datasets.narma30(n_samples=4, sample_len=2*neurons)

    # construct individual nodes
    reservoir = delay_nodes.OpticDelayNode(output_dim=neurons, beta=0.683, gamma=0.073, w=w)
    #reservoir = Oger.nodes.ReservoirNode(inputs, 100, input_scaling=0.05)
    readout = Oger.nodes.RidgeRegressionNode(ridge_param=0.000000000000000005)

    # build network with MDP framework
    flow = Oger.nodes.InspectableFlow([reservoir, readout], verbose=1)
    
    data = [x[0:-1], zip(x[0:-1], y[0:-1])]
    
    # train the flow 
    flow.train(data)
    
    #apply the trained flow to the training data and test data
    trainout = flow(x[0])
    testout = flow(x[3])

    print "NRMSE: " + str(Oger.utils.nrmse(y[3], testout))

    #plot the input
    pylab.subplot(nx, ny, 1)
    pylab.plot(x[0])
    
    #plot everything
    pylab.subplot(nx, ny, 2)
    pylab.plot(trainout, 'r')
    pylab.plot(y[0], 'b')

    pylab.subplot(nx, ny, 3)
    pylab.plot(testout, 'r')
    pylab.plot(y[3], 'b')
    
    pylab.subplot(nx, ny, 4)
    pylab.plot(flow.inspect(reservoir))
    pylab.show()
    
