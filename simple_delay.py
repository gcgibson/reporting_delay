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
from pydlm import dlm, trend, seasonality, dynamic, autoReg, longSeason

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
        if len(elm[1]) == 1:
            elm[1] = "0" + elm[1]
        date_to_index[elm[0]+elm[1]] = i
    
        i+=1
    
    d_to_i = {}
    i = 0
    iter_ =  date_to_index.keys()
    iter_.sort()
    for key in iter_:
    
        d_to_i[key] = i
        i+=1
    
    n_t = np.zeros((52-1,52-1))
    
    for elm in n_t_d:
        try:
            
            sick_date = d_to_i[elm[0]+elm[1]]
            report_date = d_to_i[elm[4] + elm[5]]
            if int(elm[4] + elm[5]) < 201621 and int(elm[3]) == 1:
                n_t[sick_date][report_date] += int(elm[3])
            
        except:
            pass
    
    
    D = 4
    
    n_t_d = []
    for row in range(len(n_t)):
        #if len(n_t[row][row:row+D]) == D:
            tmp = n_t[row][row:row+D].tolist()
            while len(tmp) < D:
                tmp += [0]
            n_t_d.append(tmp)
    
    n_t_d = np.array(n_t_d)
    
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
        
    

train_n_t_d = n_t_d[:len(n_t_d)-D + 1]
train_n_t_inf = n_t_inf[:len(n_t_d)-D + 1]

test_n_t_d = n_t_d[len(n_t_d)- D + 1:]
test_n_t_inf = n_t_inf[len(n_t_d)-D +1:]

import numpy as np
import pymc3 as pm
import pandas as pd


## Delay Model

DELAY_DIST = True
if DELAY_DIST == True:

    k = np.array(train_n_t_d).shape[1 ]
    
    with pm.Model() as multinom_test:
        a = pm.Dirichlet('a', a=np.ones(k))
        for i in range(len(train_n_t_d)):
            data_pred = pm.Multinomial('data_pred_%s'% i, n=sum(train_n_t_d[i]), p=a, observed=train_n_t_d[i])
        trace = pm.sample(50000, pm.Metropolis())
        #trace = pm.sample(1000) # also works with NUTS
    
    pm.traceplot(trace[500:]);

state_trajectories = []
PF = False
if  PF:
    N = 10000
    state_space_dimension = 1
    D = 4
    params = []
    means , particles, weights = run_pf(train_n_t_inf,N,state_space_dimension,D,params)
    
    
    
    ### Interval Predictions
    state_trajectories = [particles]
    observation_trajectories = [np.exp(particles)]
    for i in range(len(test_n_t_inf)):
        tmp = expected_value_transition_function(state_trajectories[i-1])
        observation_trajectories.append(expected_value_observation_function(tmp))
        state_trajectories.append(tmp) 
    
    state_trajectories = state_trajectories[1:]
    ## MEAN
    print (np.mean(observation_trajectories,axis=1))
    ## QUANTILES 
    state_trajectories = np.array(state_trajectories).reshape((len(test_n_t_inf),-1))


else:
    myDLM = dlm(train_n_t_inf)
    myDLM = myDLM + trend(1, name='lineTrend', w=1.0)
    # add a 7 day seasonality with prior covariance 1.0
    myDLM = myDLM + seasonality(52, name='7day', w=1.0)
    # add a 3 step auto regression
    myDLM = myDLM + autoReg(degree=2, data=train_n_t_inf, name='ar3', w=1.0)
    myDLM.fit()
    (predictMean, predictVar) = myDLM.predictN(N=3, date=myDLM.n-1)

     


for i in range(len(predictMean)):
    samples = np.random.normal(predictMean[i],np.sqrt(predictVar[i]),4)
    state_trajectories.append(samples)
state_trajectories = np.array(state_trajectories)


phat = trace['a'].mean(axis=0)
from scipy.stats import binom



myDLM.plot()

##compute weighted trajectories 

weighted_trajectories = []
for i in range(len(state_trajectories)):
    tmp = []
    samples = state_trajectories[i]
    row_sum = sum(test_n_t_d[i])
    q = sum(phat[:len(phat)-i-1])
    for samp in samples:
        btemp = binom.pmf(row_sum,samp,q)
        if np.isnan(btemp):
            tmp.append(0)
        else:
            tmp.append(btemp)
        print (row_sum,samp,q,btemp)
    weighted_trajectories.append(tmp)
weighted_trajectories = np.array(weighted_trajectories)


for i in range(len(weighted_trajectories)):
    weighted_trajectories[i] = weighted_trajectories[i]/sum(weighted_trajectories[i])


###
from sklearn.metrics import mean_squared_error

print (mean_squared_error(np.average(state_trajectories,axis=1),test_n_t_inf))
print (mean_squared_error(np.average(state_trajectories,weights = weighted_trajectories,axis=1),test_n_t_inf))


max_indeces = np.argmax(weighted_trajectories,axis=1)
max_point = []
for i in range(len(max_indeces)):
    max_point.append(state_trajectories[i][max_indeces[i]])



print (mean_squared_error(max_point,test_n_t_inf))






