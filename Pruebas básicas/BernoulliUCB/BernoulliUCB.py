# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 14:32:25 2022

@author: marco
"""

import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

def cumRegretUCB(n,k,delta,arms,gaps,nsamples):
    cumRegret = np.zeros(n)
    
    for i in range(nsamples):
        T = k*[0] # número de veces que se ha elegido cada brazo
        meanReward = k*[0] # media muestral de las recompensas obtenidas por cada brazo
        uCB = np.full(k,np.inf) # cota superior de confianza de cada brazo
        regret = 0
      
        for t in range(n):
            
            chosenArm = max(range(k),key=lambda i: uCB[i])
            rwd = arms[chosenArm].rvs()
            meanReward[chosenArm] = T[chosenArm]/(T[chosenArm]+1)*meanReward[chosenArm] \
                                    + rwd/(T[chosenArm]+1)
            T[chosenArm] +=1
            
            uCB[chosenArm] = meanReward[chosenArm] + math.sqrt((2*math.log(1/delta))/T[chosenArm])
            
            regret += gaps[chosenArm]
            cumRegret[t] += regret
        
    
    cumRegret/= nsamples
    return cumRegret

n = 1000
k = 2
probs = [0.5,0.49]
arms = [stats.bernoulli(p) for p in probs]
gaps = [0,0.01]
nsamples = 1000

fig = plt.figure()

delta = 1
label = 'δ = {}'.format(delta)
cumRegret = cumRegretUCB(n,k,delta,arms,gaps,nsamples)
plt.plot(cumRegret,color='tab:purple',label=label)

delta = 0.8
label = 'δ = {}'.format(delta)
cumRegret = cumRegretUCB(n,k,delta,arms,gaps,nsamples)
plt.plot(cumRegret,color='tab:blue',label=label)

delta = 0.5
label = 'δ = {}'.format(delta)
cumRegret = cumRegretUCB(n,k,delta,arms,gaps,nsamples)
plt.plot(cumRegret,color='tab:green',label=label)

delta = 0.1
label = 'δ = {}'.format(delta)
cumRegret = cumRegretUCB(n,k,delta,arms,gaps,nsamples)
plt.plot(cumRegret,color='tab:cyan',label=label)

delta = 0.01
label = 'δ = {}'.format(delta)
cumRegret = cumRegretUCB(n,k,delta,arms,gaps,nsamples)
plt.plot(cumRegret,color='tab:olive',label=label)

delta = 0.001
label = 'δ = {}'.format(delta)
cumRegret = cumRegretUCB(n,k,delta,arms,gaps,nsamples)
plt.plot(cumRegret,color='tab:red',label=label)


plt.xlabel('Ronda')
plt.ylabel('Remordimiento acumulado')
plt.legend(loc='upper left')
fig.savefig('Bernoulli3UCB1000.pdf',format='pdf')
plt.show()

