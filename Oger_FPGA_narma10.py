import Oger
import mdp
import pylab
import delay_nodes
import numpy as np

if __name__ == "__main__":

    inputs = 1
    neurons = 150
    n_samples = 4
    n_wash = 0
    
    # parameters of model
    beta = 0.683
    gamma = 0.0273
    phi0 = np.pi/4. 
    w = np.array(np.loadtxt('w.dat'), ndmin=2)
    
    [inputs, outputs] = Oger.datasets.narma10(n_samples=n_samples+n_wash, sample_len=2*neurons)
    wash=inputs[:n_wash]
    data = [inputs[n_wash:], zip(inputs[n_wash:], outputs[n_wash:])]

    # construct individual nodes
    reservoir = delay_nodes.OpticDelayNode(output_dim=neurons, beta=beta, gamma=gamma, phi0=phi0, w=w)
    readout = Oger.nodes.RidgeRegressionNode(use_pinv=False, ridge_param=0.0000000000000001)
    #readout = mdp.nodes.LinearRegressionNode(use_pinv=True)
    #readout = mdp.nodes.RidgeScikitsLearnNode() 
    
    # Wash out
    for xwash in wash:
        dump=reservoir(xwash)

    # build network with MDP framework
    #flow = Oger.nodes.InspectableFlow([reservoir, readout], verbose=1)
    flow = reservoir + readout
       
    print "2-fold cross-validation"
    print "Cross_validate_function = crossvalidation.cross_validate"
    errors = Oger.evaluation.validate(data, flow, Oger.utils.nrmse, n_folds=int(n_samples/2))
    print "NRMSE = " + str(errors)
    print "============================================\n"

    #~ print "Leave-one-out cross-validation"
    #~ errors = Oger.evaluation.validate(data, flow, Oger.utils.nrmse, cross_validate_function=Oger.evaluation.leave_one_out)
    #~ print "NRMSE = " + str(errors)
    
    
    # train the flow 
    #flow.train(data)
    
    #~ #apply the trained flow to the training data and test data
    #~ trainout = flow(x[0])
    #~ testout = flow(x[3])
#~ 
    #~ print "NRMSE: " + str(Oger.utils.nrmse(y[3], testout))
#~ 
    #~ #plot the input
    #~ pylab.subplot(3, 1, 1)
    #~ pylab.plot(x[0])
    #~ 
    #~ pylab.subplot(3, 1, 2)
    #~ pylab.plot(trainout, 'r')
    #~ pylab.plot(y[0], 'b')
#~ 
    #~ pylab.subplot(3, 1, 3)
    #~ pylab.plot(testout, 'r')
    #~ pylab.plot(y[3], 'b')
    #~ 
    #~ pylab.subplot(nx, ny, 4)
    #~ pylab.plot(flow.inspect(reservoir))
    #~ pylab.show()
    #~ 
