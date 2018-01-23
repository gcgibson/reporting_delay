#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  6 17:12:30 2018

@author: gcgibson
"""

import scipy.stats
import sys
import numpy as np
from numpy.random import random
import math

def observation_function(time_series_at_t,t,D,particle,params):
    
    tmp =  scipy.stats.poisson.pmf(time_series_at_t,np.exp(particle))
    if math.isnan(tmp):
        tmp = 0
    return tmp

def transition_function(particles,params):
    particles[:,0]  += np.random.normal(0,1,len(particles))
    return particles

def expected_value_observation_function(p):
    return np.exp(p)

def expected_value_transition_function(p):
    return p

def create_uniform_particles( N,state_space_dimension):
    particles  = np.random.normal(0,1 , size=(N,state_space_dimension))
    return particles

def predict(particles,t,params):
    particles = transition_function(particles,params)
    return particles


def update(particles, weights,ts,t,D,params):
    weights.fill(1.)
    for p in range(len(particles)):
        weights[p] *= observation_function(ts[t],t,D,particles[p],params)
    weights += 1.e-300
    return weights/sum(weights)  



def neff(weights):
    return 1. / np.sum(np.square(weights))


def estimate(particles, weights):
    """returns mean and variance of the weighted particles"""
    pos = particles[:, 0]
    mean = np.average(pos, weights=weights, axis=0)
    var  = np.average((pos - mean)**2, weights=weights, axis=0)
    return mean, var

### VARIOUS RESAMPLING SCHEMES
def multinomal_resample(weights):
    cumulative_sum = np.cumsum(weights)
    cumulative_sum[-1] = 1.  # avoid round-off errors
    return np.searchsorted(cumulative_sum, random(len(weights)))


def resample_from_index(particles, weights, indexes):
    particles[:] = particles[indexes]
    weights[:] = weights[indexes]
    weights.fill(1.0 / len(weights))
    return particles,weights

def stratified_resample(weights):
    N = len(weights)
    # make N subdivisions, chose a random position within each one
    positions = (random(N) + range(N)) / N

    indexes = np.zeros(N, 'i')
    cumulative_sum = np.cumsum(weights)

    i, j = 0, 0
    while i < N:
        if positions[i] < cumulative_sum[j]:
            indexes[i] = j
            i += 1
        else:
            j += 1
    return indexes



def run_pf(time_series,N,state_space_dimension,D,params):
    
    particles = create_uniform_particles(N=N,state_space_dimension=state_space_dimension)
    weights = np.zeros(N)    
    xs = [] 
    ws = []
    ws.append(weights)
    for t in range(len(time_series)):
        particles = predict(particles,t,params)       
        # incorporate measurements
        weights = update(particles, weights,time_series, t, D, params)
        ws.append(weights)
        #print (neff(weights),time_series[t],params)
        indexes = stratified_resample(weights)
        particles,weights = resample_from_index(particles, weights, indexes)
        mu, var = estimate(particles, weights)
        xs.append(mu.tolist())
    return xs,particles,ws

 


"""
Suppose we are given data in the form 

n_{0,0} , n_{0,1} , n_{0,2}, n_{0,3}
n_{1,0} , n_{1,1} , n_{1,2}, 0
n_{2,0} , n_{2,1} , 0      , 0
n_{3,0} , 0       , 0      , 0

where T=3.

Now suppose D=2, then we truncate this matrix as follows

n_{0,0} , n_{0,1} , n_{0,2}, 0
n_{1,0} , n_{1,1} , n_{1,2}, 0
n_{2,0} , n_{2,1} , 0      , 0
n_{3,0} , 0       , 0      , 0


so for a setting of parameters (T=2,D=3) the reporting trapezoid is completely 
defined

In order to get the N_{t,T}s we simply add up the rows


"""
SIMULATE = False

if SIMULATE == False:
    n_t_d = []
    with open("province-biweek_with_delays.csv") as f:
        i = 0
        for line in f.readlines():
            if i > 0:
                n_t_d.append(line.replace("\n","").split(','))
            i+=1
    date_to_index = {}
    
    i = 0
    for elm in n_t_d:
        date_to_index[elm[0]+elm[1]] = i
    
        i+=1
    
    
    d_to_i = {}
    i = 0
    iter_ =  date_to_index.keys()
    iter_.sort()
    for key in iter_:
    
        d_to_i[key] = i
        i+=1
    
    n_t = np.zeros((i,i))
    
    for elm in n_t_d:
        try:
            sick_date = d_to_i[elm[0]+elm[1]]
            report_date = d_to_i[elm[4] + elm[5]]
            n_t[sick_date][report_date] += int(elm[3])
        except:
            pass
    
    
    D = 4
    
    n_t_d = []
    for row in range(len(n_t)):
        if len(n_t[row][row:row+D]) == D:
            n_t_d.append(n_t[row][row:row+D].tolist())
    
    n_t_d = np.array(n_t_d)
    n_t_d = n_t_d[:,0:D]
    
    n_t_inf = np.sum(n_t_d,axis=1)

else:
    true_p = [.1,.3,.2,.15,.15,.1]
    n_t_d = []
    n_t_inf = []
    states = []
    tmp = np.random.normal(4,1 ,1)
    states.append(tmp)
    n_t_inf.append(np.random.poisson(np.exp(tmp))[0])
    for i in range(1,10):
        tmp = np.random.normal(states[i-1],1,1)
        states.append(tmp)
        n_t_inf.append(np.random.poisson(np.exp(tmp))[0])
    for i in range(len(n_t_inf)):
        n_t_d.append(np.random.multinomial(n_t_inf[i],true_p).tolist())
        
    

train_n_t_d = n_t_d[:60]
train_n_t_inf = n_t_inf[:60]

test_n_t_d = n_t_d[60:]
test_n_t_inf = n_t_inf[60:]

import numpy as np
import pymc3 as pm
import pandas as pd


## Delay Model

DELAY_DIST = False
if DELAY_DIST == True:

    k = np.array(train_n_t_d).shape[1 ]
    
    with pm.Model() as multinom_test:
        a = pm.Dirichlet('a', a=np.ones(k))
        for i in range(len(train_n_t_d)):
            data_pred = pm.Multinomial('data_pred_%s'% i, n=sum(train_n_t_d[i]), p=a, observed=train_n_t_d[i])
        trace = pm.sample(50000, pm.Metropolis())
        #trace = pm.sample(1000) # also works with NUTS
    
    pm.traceplot(trace[500:]);


### Process Model
if False:
	N = 10000
	state_space_dimension = 1
	D = 6
	params = []
	means , particles, weights = run_pf(train_n_t_inf,N,state_space_dimension,D,params)

	import matplotlib.pyplot as plt



	with pm.Model() as model:
	    l = pm.Gamma("l", alpha=2, beta=1)
	    nu = pm.HalfCauchy("nu", beta=5)

	    cov = nu**2 * pm.gp.cov.Matern52(1, l)
	    gp = pm.gp.Latent(cov_func=cov)

	    f = gp.prior("f", X=np.array(range(len(train_n_t_inf))).reshape((-1,1)))

	    sigma = pm.HalfCauchy("sigma", beta=5)
	    v = pm.Gamma("v", alpha=2, beta=0.1)
	    y_ = pm.NegativeBinomial("y", mu=np.exp(f), alpha=v, observed=train_n_t_inf)

	    trace = pm.sample(10000)




	n_new = len(test_n_t_inf)
	test_n_t_inf = np.array(test_n_t_inf).reshape((-1,1))

	# add the GP conditional to the model, given the new X values
	with model:
	    f_pred = gp.conditional("f_pred", test_n_t_inf)

	# Sample from the GP conditional distribution
	with model:
	    pred_samples = pm.sample_ppc(trace, vars=[f_pred], samples=1000)


	plt.plot(np.exp(pred_samples['f_pred'].mean(axis=0)))
	plt.plot(test_n_t_inf)
	plt.show()


	import pyflux as pf
	financial_crises = pd.DataFrame(n_t_inf)
	financial_crises.index = range(len(n_t_inf))
	financial_crises.columns = ["Number of banking crises"]

	model = pf.GAS(ar=4,sc=2,data=financial_crises,family=pf.Poisson())
	x = model.fit()
	gas_preds = model.predict(len(test_n_t_inf))

	print (mean_squared_error(gas_preds,test_n_t_inf))



import matplotlib.pyplot as plt
import matplotlib.cm as cmap

import numpy as np
np.random.seed(206)

import numpy as np
import matplotlib.pyplot as pl
from scipy.integrate import odeint

# Test data
n = 10
Xtest = np.linspace(1, 10, n).reshape(-1,1)

# Define the kernel function
def kernel(a, b, param):
    sqdist = np.sum(a**2,1).reshape(-1,1) + np.sum(b**2,1) - 2*np.dot(a, b.T)
    return np.exp(-.5 * (1/param) * sqdist)

#define mean function
def mean_func(x):
    ret_array = []
    for elm in x:
        # Total population, N.
        N = 1000
        # Initial number of infected and recovered individuals, I0 and R0.
        I0, R0 = 1, 0
        # Everyone else, S0, is susceptible to infection initially.
        S0 = N - I0 - R0
        # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
        beta, gamma = 3, 1 
        # A grid of time points (in days)
        t = np.linspace(0, elm, elm)
        
        # The SIR model differential equations.
        def deriv(y, t, N, beta, gamma):
            S, I, R = y
            dSdt = -beta * S * I / N
            dIdt = beta * S * I / N - gamma * I
            dRdt = gamma * I
            return dSdt, dIdt, dRdt
        
        # Initial conditions vector
        y0 = S0, I0, R0
        # Integrate the SIR equations over the time grid, t.
        ret = odeint(deriv, y0, t, args=(N, beta, gamma))
        S, I, R = ret.T
        ret_array.append(I[-1])
        
    return np.array(ret_array).reshape((-1,1))

param = 0.1
K_ss = kernel(Xtest, Xtest, param)

# Get cholesky decomposition (square root) of the
# covariance matrix
L = np.linalg.cholesky(K_ss + 1e-15*np.eye(n))
# Sample 3 sets of standard normals for our test points,
# multiply them by the square root of the covariance matrix
f_prior = np.dot(L, np.random.normal(size=(n,3)))

# Noiseless training data
Xtrain = np.array([1,2,3,4,5]).reshape(5,1)
ytrain = np.sin(Xtrain)

# Apply the kernel function to our training points
K = kernel(Xtrain, Xtrain, param)
L = np.linalg.cholesky(K + 0.00005*np.eye(len(Xtrain)))

# Compute the mean at our test points.
K_s = kernel(Xtrain, Xtest, param)
Lk = np.linalg.solve(L, K_s)
mu = mean_func(Xtest) + np.dot(Lk.T, np.linalg.solve(L, (ytrain-mean_func(Xtrain)))).reshape((n,1))

mu = mu.reshape((n,))
# Compute the standard deviation so we can plot it
s2 = np.diag(K_ss) - np.sum(Lk**2, axis=0)
stdv = np.sqrt(s2)
# Draw samples from the posterior at our test points.
L = np.linalg.cholesky(K_ss + 1e-6*np.eye(n) - np.dot(Lk.T, Lk))
f_post = mu.reshape(-1,1) + np.dot(L, np.random.normal(size=(n,3)))

pl.plot(Xtrain, ytrain, 'bs', ms=8)
pl.plot(Xtest, f_post)
pl.gca().fill_between(Xtest.flat, mu-2*stdv, mu+2*stdv, color="#dddddd")
pl.plot(Xtest, mu, 'r--', lw=2)
pl.title('Three samples from the GP posterior')
pl.show()


















sys.exit()

print (np.exp(means).tolist())
print (train_n_t_inf.tolist())
mean_squared_error(np.exp(means).tolist(),train_n_t_inf.tolist())
mean_squared_error(gpy_preds,test_n_t_inf.tolist())
### Point Predictions
p_d_hat = np.mean(trace['a'],axis=0)
n_t_inf_hat = []
n_t_d_hat = []
for i in range(len(test_n_t_inf)):
    tmp_n_t_inf = np.exp(means[len(means)-1])
    n_t_inf_hat.append(tmp_n_t_inf)
    n_t_d_hat.append((tmp_n_t_inf*p_d_hat).tolist())
    


print mean_squared_error(n_t_d_hat,test_n_t_d)


### Interval Predictions
state_trajectories = [particles]
observation_trajectories = [np.exp(particles)]
for i in range(1,len(test_n_t_inf)):
    tmp = expected_value_transition_function(state_trajectories[i-1])
    observation_trajectories.append(expected_value_observation_function(tmp))
    state_trajectories.append(tmp) 



## MEAN
print (np.mean(observation_trajectories,axis=1))
## QUANTILES 
print (np.percentile(observation_trajectories,95,axis=1))
print (np.percentile(observation_trajectories,5,axis=1))






