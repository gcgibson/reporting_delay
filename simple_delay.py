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
    
    tmp =  scipy.stats.poisson.pmf(time_series_at_t,particle)
    if math.isnan(tmp):
        tmp = 0
    return tmp
def transition_function(particles,params):
    
    particles[:,0]  += np.random.normal(0,10,len(particles))
    return particles


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
        xs.append(mu)
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
SIMULATE = True

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
        
    
    

import numpy as np
import pymc3 as pm
import pandas as pd

DELAY_DIST = False
if DELAY_DIST == True:

    k = np.array(n_t_d).shape[1 ]
    
    with pm.Model() as multinom_test:
        a = pm.Dirichlet('a', a=np.ones(k))
        for i in range(sample_size):
            data_pred = pm.Multinomial('data_pred_%s'% i, n=sum(n_t_d[i]), p=a, observed=n_t_d[i])
        trace = pm.sample(50000, pm.Metropolis())
        #trace = pm.sample(1000) # also works with NUTS
    
    pm.traceplot(trace[500:]);
N = 10000
state_space_dimension = 1
D = 6
params = []
print (n_t_inf)
means , particles,weights = run_pf(n_t_inf,N,state_space_dimension,D,params)


print (means)
print (states)



