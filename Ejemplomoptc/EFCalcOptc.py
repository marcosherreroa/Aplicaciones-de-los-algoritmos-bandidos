# -*- coding: utf-8 -*-
"""
Aplicaciones de los algoritmos de bandidos. Autor: Marcos Herrero
Bandidos estocásticos
Algoritmo Explora-Primero
Ejemplo m óptimo c)
"""

import math
import random
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np

def regretEF(n,k,m,arms,gaps):
    rwds = k*[0]
    
    for i in range(m):
        for j in range(k):
            rwds[j] += arms[j].rvs()
    
    maximum = max(rwds)
    bestarm = random.choice([i for i in range(k) if rwds[i] == maximum])
    
    return m*sum(gaps)+(n-m*k)*gaps[bestarm]

n = 200
Delta = 0.3
arms = [stats.norm(0,1),stats.norm(-Delta,1)]
gaps = [0,Delta]

expectedRegret = np.empty(n//2+1)
expectedRegret[0] = 0.5*n*Delta

X = stats.norm(0,1)

for m in range(1,n//2+1):
    expectedRegret[m] = m*Delta+(n-m)*Delta*X.cdf(-m*Delta/math.sqrt(2*m))


ms = np.arange(0,n//2+1,10)
nsamples = 100
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
plt.plot(expectedRegret,color='tab:blue',label='Calculado a partir de la distribución',linewidth = 4)
plt.plot(ms,averageRegret,color='tab:green',label='Estimado con {} simulaciones'.format(nsamples))
plt.xlabel('m')
plt.ylabel('Remordimiento esperado')
plt.legend(loc='upper right')
fig.savefig('EFCalcOptc.pdf',format='pdf')
plt.show()
