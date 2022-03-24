# -*- coding: utf-8 -*-
"""
Experimento ETC
Autor : Marcos Herrero
"""

import math
import random
import sympy.stats as stats
import matplotlib.pyplot as plt
import numpy as np

def computem (n,Delta):
    if Delta == 0:
        return 0
    else :
        return max(1,math.ceil(4/(Delta*Delta)*math.log(n*Delta*Delta/4)))

def regretExploreFirst(n,k,m,arms,gaps):
    rwds = k*[0]
    
    for i in range(m):
        for j in range(k):
            rwds[j] += stats.sample(arms[j])
    
    maximum = max(rwds)
    bestarm = random.choice([i for i in range(k) if rwds[i] == maximum])
    #bestarm = max(range(k),key=lambda i:rwds[i])
    
    return m*sum(gaps)+(n-m*k)*gaps[bestarm]

def plotDeltaRegret():
       
    n = 1000
    sampleNum = 1000
    
    arms = 2*[None]
    arms[0] =  stats.Normal('X0',0,1)
    gaps = 2*[0]
    
    nDeltas = 40
    Deltas = np.linspace(0,1,nDeltas)
    regret = np.empty(nDeltas)
    upperBound = np.empty(nDeltas)
    
    for i in range(nDeltas):
        Delta = Deltas[i]
        
        arms[1]= stats.Normal('X1',-Delta,1)
        m = computem(n,Delta)
        gaps[1] = Delta
        
        regret[i] = 0
        for k in range(sampleNum):
            regret[i] += regretExploreFirst(n,2,m,arms,gaps)
        
        regret[i] /= sampleNum
        
        if Delta == 0:
            upperBound[i] = 0
        else:
            upperBound[i] = min(n*Delta,Delta+4/Delta*(1+max(0,math.log(n*Delta*Delta/4))))
    
    plt.plot(Deltas,regret, color='tab:blue')
    plt.plot(Deltas,upperBound,'--', color='tab:blue')
    plt.xlabel('∆')
    plt.ylabel('Remordimiento')
    fig = plt.gcf()
    fig.savefig('EFExpDeltaRegret.pdf',format='pdf')
    plt.show()

def plotmRegretAndVar():
    
    n = 1000
    Delta = 0.1
    sampleNum = 600
    
    arms = 2*[None]
    arms[0] =  stats.Normal('X0',0,1)
    arms[1] = stats.Normal('X1',-Delta,1)
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
            regretSamples[k] = regretExploreFirst(n,2,m,arms,gaps)
        
        regret[i] = sum(regretSamples)/sampleNum
        
        regretDev[i] = 0
        for k in range(sampleNum):
            dif = regretSamples[k]-regret[i]
            regretDev[i] += dif*dif
        
        regretDev[i] = math.sqrt(regretDev[i]/(n-1))
    
    indmin = min(range(numms), key = lambda i: regret[i])
    print('Valor de m con regret mínimo: {}'.format(ms[indmin]))
    print('Regret mínimo: {}'.format(regret[indmin]))
      
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
    
#plotDeltaRegret()
plotmRegretAndVar()