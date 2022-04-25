# -*- coding: utf-8 -*-
"""
Aplicaciones de los algoritmos de bandidos. Autor: Marcos Herrero
Concentración de la medida
Comparación del tamaño de las colas con la cota proporcionada por Berry-Esseen-Mills
"""

import math
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np

#Compara media muestral de Bernoullis
def comparaBernoulli():
    #Props Bernoulli
    p = 0.5
    mu = p
    sigma = math.sqrt(p*(1-p))
    rho = p*(1-p)*(p**2+(1-p)**2)
    
    numns = 1000
    epsilon = 0.25*sigma
    
    chebBound = np.empty(numns)
    CLTBound = np.empty(numns)
    realTail = np.empty(numns)
    
    for i in range(numns):
        n = i+1
        
        chebBound[i] = sigma**2/(n*(epsilon**2))
        CLTBound[i] = math.sqrt(2/math.pi)*math.exp(-n*epsilon**2/(2*sigma**2)) + \
            rho/(sigma**3*math.sqrt(n))
            
        B = stats.binom(n,p)
        realTail[i] = 1 - B.cdf(math.ceil(n*(mu+epsilon)-1)) + B.cdf(n*(mu-epsilon))
    
    return np.array(range(1,numns+1)),chebBound,CLTBound,realTail

def comparaPoisson():
    lamb = 100
    mu = lamb
    sigma = math.sqrt(lamb)
    rho = lamb
    
    numns = 1000
    epsilon = 0.25*sigma
    
    chebBound = np.empty(numns)
    CLTBound = np.empty(numns)
    realTail = np.empty(numns)
    
    for i in range(numns):
        n = i+1
        
        chebBound[i] = sigma**2/(n*(epsilon**2))
        CLTBound[i] = math.sqrt(2/math.pi)*math.exp(-n*epsilon**2/(2*sigma**2)) + \
            rho/(sigma**3*math.sqrt(n))
            
        P = stats.poisson(n*lamb)
        realTail[i] = 1 - P.cdf(math.ceil(n*(mu+epsilon)-1)) + P.cdf(n*(mu-epsilon))
    
    return np.array(range(1,numns+1)),chebBound,CLTBound,realTail

ns,chebBound,CLTBound, realTail = comparaBernoulli()


fig = plt.figure()
ax = plt.subplot(111)
plt.plot(ns,chebBound, '--',color = 'tab:green',label='Cota de Chebyshev')
plt.plot(ns,CLTBound, '--',color = 'tab:blue',label='Cota Berry-Esseen + Mills')
plt.plot(ns,realTail, color = 'tab:blue',label='P(|X-μ| ≥ ε')
plt.xlabel('n')
plt.ylim([0,0.05])
plt.legend(loc='upper right')
fig.savefig('CompBernoulliCLT.pdf',format='pdf')
plt.show()    
                                                                   
#                                                                   
