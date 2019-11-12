from numpy import sin, sum, arange, zeros, size, ones, ceil, \
                  histogram, round, hstack, transpose, dot, \
                  std, sqrt 
from numpy.random import rand 
from scipy import exp, pi 


def optic(x,xd,al,bet):
    phi0=0.25*pi
    return al*x+bet*(sin(xd+phi0)**2)

def narma10(nar):
    ord=10
    warm=ord*30
    nar=nar+warm
    
    # order 10
    a=.3
    b=.05
    c=1.5
    d=.1 # 0.1 -- 0.001
    
    inputs = rand(nar,)*0.5
    outputs = 0.1*ones(nar)
    for i in range(ord,nar):
        outputs[i] = a*outputs[i-1] + b*outputs[i-1]*sum(outputs[i-ord:i-1]) \
        + c*inputs[i-ord] * inputs[i-1] + d
        
    inputs=inputs[warm:]
    outputs=outputs[warm:]
    return inputs, outputs

def narma10map(nar,b):
    ord=10
    warm=ord
    nar=nar+warm
    
    # order 10
    a=.3
    #b=.05
    #c=1.5
    d=.1 # 0.1 -- 0.001
    
    #inputs = rand(nar,)*0.5
    outputs = 0.1*ones(nar)+rand(nar,)*0.2
    for i in range(ord,nar):
        outputs[i] = a*outputs[i-1] + b*outputs[i-1]*sum(outputs[i-ord:i-1]) \
        + d
        
    #inputs=inputs[warm:]
    outputs=outputs[warm:]
    return outputs
