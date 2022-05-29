# -*- coding: utf-8 -*-

'''
Bandidos estocásticos: introducción, algoritmos y experimentos
TFG Informática
Sección 7.2.2
Figura 6
Autor: Marcos Herrero Agustín
'''

import math
import random
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np

def computemExpl(n):
    return math.ceil(n**(2/3))

def computemTeor (n,Delta):
    if Delta == 0:
        return 0
    else :
        return max(1,math.ceil(4/(Delta*Delta)*math.log(n*Delta*Delta/4)))
    
def computemOpt(n,Delta):
    expectedRegret = np.empty(n//2+1)
    X = stats.norm(0,1)
    
    expectedRegret[0] = 0.5*n*Delta
    for m in range(1,n//2+1):
        expectedRegret[m] = m*Delta+(n-m)*Delta*X.cdf(-m*Delta/math.sqrt(2*m))
    
    mOpt = min(range(n//2+1),key = lambda i: expectedRegret[i])
    return mOpt

def samplePseudoRegretEF(n,k,m,arms,gaps):
    rwds = k*[0]
    
    for i in range(m):
        for j in range(k):
            rwds[j] += arms[j].rvs()
    
    maximum = max(rwds)
    bestarm = random.choice([i for i in range(k) if rwds[i] == maximum])
    
    return m*sum(gaps)+(n-m*k)*gaps[bestarm]

n = 1000
Delta = 0.1
sampleNum = 600

arms = 2*[None]
arms[0] = stats.norm(0,1)
arms[1] = stats.norm(-Delta,1)
gaps = 2*[0]
gaps[1] = Delta

ms = np.arange(0,401,40)
numms = ms.size
regret = np.empty(numms)
regretDev = np.empty(numms)


for i in range(numms):
    
    m = ms[i]
    
    regretSamples = np.empty(sampleNum)
    for k in range(sampleNum):
        regretSamples[k] = samplePseudoRegretEF(n,2,m,arms,gaps)
    
    regret[i] = sum(regretSamples)/sampleNum
    
    regretDev[i] = 0
    for k in range(sampleNum):
        dif = regretSamples[k]-regret[i]
        regretDev[i] += dif*dif
    
    regretDev[i] = math.sqrt(regretDev[i]/(sampleNum-1))

indmin = min(range(numms), key = lambda i: regret[i])
print('Valor de m con regret mínimo: {}'.format(ms[indmin]))
print('Regret esperado mínimo: {}'.format(regret[indmin]))
print('Varianza: {}'.format(regretDev[indmin]))
  
plt.plot(ms,regret, color='tab:blue')
plt.xlabel('m')
plt.ylabel('Remordimiento')
fig = plt.gcf()
fig.savefig('EFExpmRegret.pdf',format='pdf')
plt.show()

plt.plot(ms,regretDev, color='tab:blue')
plt.xlabel('m')
plt.ylabel('Desviación típica del remordimiento')
fig = plt.gcf()
fig.savefig('EFExpmRegretDev.pdf',format='pdf')
plt.show()
