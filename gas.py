#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 07:53:33 2018

@author: gcgibson
"""

"""
Generate GAS simulated data
"""
import numpy as np
from scipy.stats import poisson
import sys
A = .101
B = -.395
w = .183

t = 0
f = [0]
y = [np.random.poisson(np.exp(f[0]),1)[0]]
size = 100

s = []
for t in range(1,size):
    s_t_1 = (y[t-1]/np.exp(f[t-1])-1)
    s.append(s_t_1)
    f_t = w + A*s_t_1 + B* f[t-1]
    f.append(f_t)
    y.append(np.random.poisson(np.exp(f[t]),1)[0])
    
    
"""
Estimate GAS Model via ML
"""


def gas_poisson(q,y=y):
    A = q[0]
    B = q[1]
    w = q[2]
    
    n = len(y)
    s = []
    f = []
    
    f.append(0)
    
    for t in range(1,n):
        s_t_1 = (y[t-1]/np.exp(f[t-1])-1)
        s.append(s_t_1)
        f_t  = w + A*s_t_1 + B*f[t-1]
        f.append(f_t)


    """
    Log-likelihood
    """
    ll_sum = 0
    for t in range(n):
        ll_sum += poisson.logpmf(y[t],np.exp(f[t]))
        
    return (-1.0*ll_sum/n )


from scipy import optimize
#result = optimize.minimize(gas_poisson,[1,1,1],args=(y,),method="BFGS")
result = optimize.minimize(gas_poisson,[.001,.001,.001],args=(y,),method="",tol=1e-10)
