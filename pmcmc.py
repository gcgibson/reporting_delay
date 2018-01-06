from numpy.linalg import norm
from numpy.random import randn
import scipy.stats
import sys
import numpy as np
from numpy.random import random
from sklearn import preprocessing
      

def observation_function(time_series_at_t,particle):
    return scipy.stats.norm.pdf(time_series_at_t,particle,1)

def transition_function(particles,params):
    
    return particles  + np.random.normal(0,params[0],len(particles)).reshape((len(particles),1))

def forecast_n_ahead(estimated_states,n_ahead):
    ##
    last_state_mean = estimated_states[-1]

    last_state_mean_lag_1 = estimated_states[-2]
    forecasts_states = []
    forecasts_obs = []
    for i in range(n_ahead):
        forecasts_states.append(phi1*last_state_mean + phi2*last_state_mean_lag_1)  
        forecasts_obs.append(np.exp(forecasts_states[i]))
        last_state_mean_lag_1 = last_state_mean
        last_state_mean = forecasts_states[i]
    return forecasts_obs


def create_uniform_particles( N,state_space_dimension):
    particles  = np.random.uniform(1, 100, size=(N,state_space_dimension))
    return particles


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


def predict(particles,t,params):
    particles = transition_function(particles,params)
    return particles

def neff(weights):
    return 1. / np.sum(np.square(weights))


def update(particles, weights,ts,t):
    weights.fill(1.)
    for p in range(len(particles)):
        weights[p] *= observation_function(ts[t],particles[p][0])
    weights += 1.e-300
    return weights/sum(weights)  


def estimate(particles, weights):
    """returns mean and variance of the weighted particles"""

    pos = particles[:, 0]
    mean = np.average(pos, weights=weights, axis=0)
    var  = np.average((pos - mean)**2, weights=weights, axis=0)
    return mean, var


def multinomal_resample(weights):
    cumulative_sum = np.cumsum(weights)
    cumulative_sum[-1] = 1.  # avoid round-off errors
    return np.searchsorted(cumulative_sum, random(len(weights)))


def resample_from_index(particles, weights, indexes):
    particles[:] = particles[indexes]
    weights[:] = weights[indexes]
    weights.fill(1.0 / len(weights))
    return particles,weights

def run_pf(time_series,N,state_space_dimension,params):
    
    particles = create_uniform_particles(N=N,state_space_dimension=state_space_dimension)
    weights = np.zeros(N)    
    xs = [] 
    ws = []
    ws.append(weights)
    for t in range(len(time_series)):
        particles = predict(particles,t,params)       
        # incorporate measurements
        weights = update(particles, weights,time_series, t)
        ws.append(weights)
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
    
    
    
def particle_mcmc(time_series,num_iters,state_space_dimension):
    theta_current = 1.
    proposal_width = 10.
    posterior = [theta_current]
    for i in range(num_iters):
        # suggest new position
        theta_proposal = scipy.stats.norm(theta_current, proposal_width).rvs()
        
        params=[np.exp(theta_current)]
        # Compute likelihood by multiplying probabilities of each data point
        estimated_states, particles, ws = run_pf(time_series,N=100,state_space_dimension=state_space_dimension,params=params)   
        likelihood_current = get_log_likelihood_from_particle_filter(ws, params)
        
        params=[np.exp(theta_proposal)]
        estimated_states, particles, ws = run_pf(time_series,N=100,state_space_dimension=state_space_dimension,params=params)   
        likelihood_proposal = get_log_likelihood_from_particle_filter(ws, params)
        
        # Compute prior probability of current and proposed mu        
        prior_current = np.log(scipy.stats.norm(0, 10).pdf(theta_current))
        prior_proposal = np.log(scipy.stats.norm(0, 10).pdf(theta_proposal))
        
        p_current = likelihood_current + prior_current
        p_proposal = likelihood_proposal +  prior_proposal
        
        # Accept proposal?
        p_accept = np.exp(p_proposal -  p_current)
        
        # Usually would include prior probability, which we neglect here for simplicity
        accept = np.random.rand() < p_accept
        
        
        if accept:
            # Update position
            theta_current = theta_proposal
        
        posterior.append(theta_current)
        
    return posterior


time_series  = np.ones(10)
num_iters= 1000
state_space_dimension = 1

posterior = particle_mcmc(time_series,num_iters,state_space_dimension)
posterior = np.array(posterior[100:])
import matplotlib.pyplot as plt

#print
print (posterior)
print (posterior.mean())