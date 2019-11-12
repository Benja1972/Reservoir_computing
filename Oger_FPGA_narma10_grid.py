import Oger
import mdp
import pylab
import delay_nodes
import numpy as np

if __name__ == "__main__":
    
    # TODO: check influance of initialization
    
    neurons = 150
    n_samples = 4
    n_wash = 0
    
    # parameters of model
    beta = mdp.numx.arange(0.2, 1.2, 0.3)
    gamma = mdp.numx.arange(0.05, 0.55, 0.2)
    phi0 = np.pi/4. 
    w = np.array(np.loadtxt('w.dat'), ndmin=2)
    
    [inputs, outputs] = Oger.datasets.narma10(n_samples=n_samples+n_wash, sample_len=2*neurons)
    wash=inputs[:n_wash]
    data = [inputs[n_wash:], zip(inputs[n_wash:], outputs[n_wash:])]

    # construct individual nodes
    reservoir = delay_nodes.OpticDelayNode(output_dim=neurons, beta=beta[0], gamma=gamma[0], phi0=phi0, w=w)
    readout = Oger.nodes.RidgeRegressionNode(use_pinv=False, ridge_param=0.000000000001)
    #readout = mdp.nodes.LinearRegressionNode(use_pinv=True)

    # Wash out
    for xwash in wash:
        dump=reservoir(xwash)

    # build network with MDP framework
    #flow = Oger.nodes.InspectableFlow([reservoir, readout], verbose=1)
    flow = reservoir + readout
    
    mdp.caching.set_cachedir("/tmp/my_cache")
    mdp.activate_extension("cache_execute")
    
    # 2D plotting example
    print "Then we range over both beta and gamma, instantiating 2 reservoirs for each setting."
    gridsearch_parameters = {reservoir:{'beta': beta,
                                        'gamma': gamma,
                                        '_instance':range(12)}}
    print "gridsearch_parameters = " + str(gridsearch_parameters)
    opt2D = Oger.evaluation.Optimizer(gridsearch_parameters, Oger.utils.nrmse)

    # Run the gridsearch
    errors = opt2D.grid_search(data, flow, n_folds=int(n_samples/2), cross_validate_function=Oger.evaluation.n_fold_random)
    
    # Plot the results, after taking the mean over the reservoir instances
    opt2D.plot_results([(reservoir, '_instance')])
       
    #print "2-fold cross-validation"
    #print "Cross_validate_function = crossvalidation.cross_validate"
    #errors = Oger.evaluation.validate(data, flow, Oger.utils.nrmse, n_folds=int(n_samples/2))
    #print "NRMSE = " + str(errors)
    #print "============================================\n"


