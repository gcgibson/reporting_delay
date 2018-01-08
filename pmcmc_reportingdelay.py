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
    
    tmp =  scipy.stats.multinomial.pmf(time_series_at_t,particle,params)
    if math.isnan(tmp):
        tmp = 0
    return tmp
def transition_function(particles,params):
    
    particles[:,0]  += np.random.normal(0,1,len(particles))
    particles[:,1]  = np.random.poisson(np.exp(particles[:,0]))
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
        weights[p] *= observation_function(ts[t],t,D,particles[p][1],params)
    weights += 1.e-300
    return weights/sum(weights)  



def neff(weights):
    return 1. / np.sum(np.square(weights))


def estimate(particles, weights):
    """returns mean and variance of the weighted particles"""

    pos = particles[:, 1]
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

 

def get_log_likelihood_from_particle_filter(ws,params):
    log_lik = 0
    for i in range(len(ws)):
        tmp = np.array(ws[i])
        log_lik += np.log(tmp.sum()/len(tmp))
        
    return log_lik
    
    
    
def particle_mcmc(time_series,num_iters,state_space_dimension,D):
    
    
    theta_current = np.array([1./3.,1./3.,1./3.])
    posterior = [theta_current]
    
    acceptance_ratio = 0
    for i in range(num_iters):
        # suggest new position
        theta_proposal = np.random.dirichlet(100*theta_current)
        
        # Compute likelihood by multiplying probabilities of each data point
        estimated_states, particles, ws = run_pf(time_series,N=500,state_space_dimension=state_space_dimension,D=D,params=[theta_current])   
        likelihood_current = get_log_likelihood_from_particle_filter(ws, params)
        
        estimated_states, particles, ws = run_pf(time_series,N=500,state_space_dimension=state_space_dimension,D=D,params=[theta_proposal])   
        likelihood_proposal = get_log_likelihood_from_particle_filter(ws, params)
        
        
        # Compute prior probability of current and proposed mu        
        prior_current = 0# np.log(1./theta_current.sum())#np.log(scipy.stats.norm(0, 10).pdf(theta_current))
        prior_proposal =0#np.log( 1/theta_proposal.sum())#np.log(scipy.stats.norm(0, 10).pdf(theta_proposal))
        
        p_current = likelihood_current + prior_current
        p_proposal = likelihood_proposal +  prior_proposal
        
        # Accept proposal?
        p_accept = np.exp(p_proposal -  p_current)
        
        rand_ = np.random.rand()
        # Usually would include prior probability, which we neglect here for simplicity
        accept = rand_ < p_accept
        
        #print (likelihood_current,likelihood_proposal,p_accept,rand_)

        if accept:
            # Update position
            acceptance_ratio +=1
            theta_current = theta_proposal
        
        posterior.append(theta_current)
    print ("Acceptance Ratio: " + str((1.0*acceptance_ratio)/num_iters))
    return posterior



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
true_p_ds = [.25,.5,.25]

n_t_d = np.random.multinomial(10,true_p_ds,100)

#n_t_d[len(n_t_d)-1][1] = 0 
#n_t_d[len(n_t_d)-1][2] = 0
#n_t_d[len(n_t_d)-2][2] = 0 

#N_t_T = np.sum(n_t_d,axis=1)

num_iters= 100
state_space_dimension = 2
D=3


posterior = particle_mcmc(n_t_d,num_iters,state_space_dimension,D)
posterior = np.array(posterior[10:])
#
#print
print (posterior.mean(axis=0))



#means,particles,weights = run_pf(  N_t_T,N=500,state_space_dimension=state_space_dimension,D=D,params=[posterior.mean(axis=0)])

















