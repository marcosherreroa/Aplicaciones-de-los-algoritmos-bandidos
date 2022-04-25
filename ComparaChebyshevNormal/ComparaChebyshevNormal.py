# -*- coding: utf-8 -*-
"""
Aplicaciones de los algoritmos de bandidos. Autor: Marcos Herrero
Concentración de la medida
Comparación Markov, Chebyshev con N(0,1)
"""

import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np

epsilons = np.arange(2.5,6,0.2)
chebBound = np.empty(epsilons.size)
realProb = np.empty(epsilons.size)

X = stats.norm(0,1)


for i  in range(epsilons.size):
    
    chebBound[i] = 1/(epsilons[i]**2)
    realProb[i] = 2*X.cdf(-epsilons[i])

fig = plt.figure()
ax = plt.subplot(111)
plt.plot(epsilons,1/chebBound, '--',color = 'tab:blue',label='ε²')
plt.plot(epsilons,1/realProb, color = 'tab:blue',label='1/P(|X| ≥ ε)')
plt.xlabel('ε')
plt.legend(loc='upper left')
fig.savefig('CompChebNorm256.pdf',format='pdf')

plt.show()