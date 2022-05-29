# -*- coding: utf-8 -*-

"""
Bandidos estocásticos: introducción, algoritmos y experimentos
TFG Informática
Sección 7.2.7
Figura 12
Autor: Marcos Herrero Agustín
"""

import math
import random
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np


class DegenerateRV:
    
    def __init__(self,value):
        self.value = value
    
    def rvs(self):
        return self.value
    
prob = 0

def regretEF(n,k,m,arms,gaps):
    rwds = k*[0]
    
    for i in range(m):
        for j in range(k):
            rwds[j] += arms[j].rvs()
    
    maximum = max(rwds)
    bestarm = random.choice([i for i in range(k) if rwds[i] == maximum])
    global prob
    prob = prob + bestarm
    
    return m*sum(gaps)+(n-m*k)*gaps[bestarm]

n = 1000
p = 0.6
q = 0.5
arms = [stats.bernoulli(p),DegenerateRV(q)]
Delta = p - q
gaps = [0,Delta]

expectedRegret = np.empty(n//2+1)


for m in range(n//2+1):
    X = stats.binom(m,p)
    expectedRegret[m] = m*Delta+(n-m)*Delta*(X.cdf(math.ceil(q*m-1)) + \
                        0.5*(X.cdf(math.floor(q*m))-X.cdf(math.ceil(q*m-1))))


ms = np.arange(0,n//2+1,10)
nsamples = 50
numms = ms.size

averageRegret = np.zeros(numms)


for i in range(numms):
    m = ms[i]
    
    for s in range(nsamples):
        averageRegret[i] = s/(s+1)*averageRegret[i]+1/(s+1)*regretEF(n,2,m,arms,gaps)
    
mopt = min(range(n//2+1),key = lambda i: expectedRegret[i])
print('m_opt = {}'.format(mopt))
print('Regret mínimo: {}'.format(expectedRegret[mopt]))   

fig = plt.figure()
ax = plt.subplot(111)
plt.plot(expectedRegret,color='tab:blue',label='Calculado a partir de la distribución', linewidth = 4)
plt.plot(ms,averageRegret,color='tab:green',label='Estimado con {} simulaciones'.format(nsamples))
plt.xlabel('m')
plt.ylabel('Remordimiento esperado')
plt.legend(loc='upper right')
fig.savefig('EFCalcOptb.pdf',format='pdf')
plt.show()



    