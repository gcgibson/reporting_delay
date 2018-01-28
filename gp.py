#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 11:12:37 2018

@author: gcgibson
"""
import numpy as np
import GPy
X = np.arange(10).reshape((-1,1))
y = np.random.poisson(5,10).reshape((-1,1))
kernel = GPy.kern.RBF(1, variance=1.0, lengthscale=1.0)

poisson_likelihood = GPy.likelihoods.Poisson()
laplace_inf = GPy.inference.latent_function_inference.Laplace

m = GPy.core.GP(X=X, Y=y, likelihood=poisson_likelihood, inference_method=laplace_inf, kernel=kernel)

m.optimize()
import matplotlib.pyplot as plt

plt.plot(m.predict( np.arange(10).reshape((-1,1)))[0])
#plt.plot(y)
plt.show()
print m