from numpy import sin, sum, arange, zeros, size, ones, ceil, \
                  histogram, round, hstack, transpose, dot, \
                  std, sqrt 
from numpy.random import rand 
from scipy import exp, pi 
from scipy.linalg import pinv, pinv2, norm
import pylab as pl
from reservoir import optic, narma10

# Physical params
betta=.03535
ap=0.5
neur=50

# RC params
Nskip=neur+500
Nwash=200
Ntrain=200
Ntest=60
Nbigsteps=Nskip+Nwash+Ntrain+Ntest

# Input data: NARMA
seeds,sprout=narma10(Nwash+Ntrain+Ntest)

ampin=0.15
inp=ampin*hstack((zeros(Nskip,),seeds))

# Connection matrix
w=rand(neur,)
w=w/sum(abs(w))
#w=0.1*ones(d,)+0.01*w1
#w=0.1*w1


x=zeros(Nbigsteps,)
x[:neur]=0.1*rand(neur,)
rc_states=zeros((Nbigsteps,neur));

for i in range(neur,Nbigsteps):
    xh=sum(w*x[i-neur:i])+inp[i]
    x[i]=optic(x[i-1],xh,ap,betta)
    rc_states[i,:]=x[i-neur+1:i+1].copy()

# Read-out, Tarin and Test
rc_train=rc_states[Nskip+Nwash+1:Nskip+Nwash+Ntrain,:].copy()
rc_test=rc_states[Nskip+Nwash+Ntrain+1:,:].copy()

tar_train=sprout[Nwash:Nwash+Ntrain-1].copy()
tar_test=sprout[Nwash+Ntrain:-1].copy()

alpha=dot(pinv2(rc_train),tar_train)

out_train=dot(rc_train,alpha)
out_test=dot(rc_test,alpha)

# Errors
NRMSE_train=norm(tar_train-out_train)/(sqrt(len(tar_train))*std(tar_train))
NRMSE_test=norm(tar_test-out_test)/(sqrt(len(tar_test))*std(tar_test))



# Figures
plot_params = {'axes.labelsize': 18,
               'text.fontsize': 20,
               'legend.fontsize': 20,
               'title.fontsize': 20,
               'xtick.labelsize': 18,
               'ytick.labelsize': 18}
          
pl.rcParams.update(plot_params)

pl.figure(10)
pl.plot(x)
pl.plot(inp)


pl.figure(11)
pl.plot(tar_train)
pl.plot(out_train)
pl.title(r'Train, $\beta=$ %1.3f, NRMSE$_{train}$= %1.4f' 
            % (betta,  NRMSE_train) )

pl.figure(12)
pl.plot(tar_test)
pl.plot(out_test)
pl.title(r'Test, $\beta=$ %1.3f, NRMSE$_{test}$= %1.4f' 
            % (betta,  NRMSE_test))


#~ pl.figure(20)
#~ pl.imshow(transpose(rc_states))
#~ pl.axis('tight')

#hbin=500
#H, xedges = histogram(x[ceil(len(x)*0.6):], bins=hbin)


#pl.figure(2)
#pl.plot(xedges[1:],H)


pl.show()
