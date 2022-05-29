# -*- coding: utf-8 -*-
"""
Sección 8
Esperanza y varianza del pseudo-remordimiento en UCB
Autor: Marcos Herrero
"""

import math
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np


def samplePseudoRegretUCB(n,k,delta,arms,gaps):#cambiar a pseudo
    T = k*[0] # número de veces que se ha elegido cada brazo
    meanReward = k*[0] # media muestral de las recompensas obtenidas por cada brazo
    UCB = k*[np.inf] # cota superior de confianza de cada brazo
    regret = 0
    
    
    for i in range(n):
        
        chosenArm = max(range(k),key=lambda i: UCB[i])
        rwd = arms[chosenArm].rvs()
        meanReward[chosenArm] = T[chosenArm]/(T[chosenArm]+1)*meanReward[chosenArm] \
                                + rwd/(T[chosenArm]+1)
        T[chosenArm] +=1
        
        UCB[chosenArm] = meanReward[chosenArm] + math.sqrt((2*math.log(1/delta))/T[chosenArm])
        
        regret += gaps[chosenArm]
    
    return regret


n = 1000
Delta = 0.1
sampleNum = 600

arms = 2*[None]
arms[0] = stats.norm(0,1)
arms[1] = stats.norm(-Delta,1)
gaps = 2*[0]
gaps[1] = Delta

xaxis = np.array(range(12))
deltas = np.array([10**(-i) for i in range(12)])
numdeltas = deltas.size
regret = np.empty(numdeltas)
regretDev = np.empty(numdeltas)


for i in range(numdeltas):
    
    delta = deltas[i]
    
    regretSamples = np.empty(sampleNum)
    for k in range(sampleNum):
        regretSamples[k] = samplePseudoRegretUCB(n,2,delta,arms,gaps)
    
    regret[i] = sum(regretSamples)/sampleNum
    
    regretDev[i] = 0
    for k in range(sampleNum):
        dif = regretSamples[k]-regret[i]
        regretDev[i] += dif*dif
    
    regretDev[i] = math.sqrt(regretDev[i]/(sampleNum-1))

indmin = min(range(numdeltas), key = lambda i: regret[i])
print('Valor de delta con regret mínimo: {}'.format(deltas[indmin]))
print('Regret esperado mínimo: {}'.format(regret[indmin]))
print('Varianza: {}'.format(regretDev[indmin]))

fig = plt.figure()
plt.plot(xaxis,regret, color='tab:blue')
plt.xlabel('-log₁₀ δ')
plt.ylabel('Remordimiento')
fig.savefig('EspPseudoregUCB.pdf',format='pdf')
plt.show()

plt.plot(xaxis,regretDev, color='tab:blue')
plt.xlabel('-log₁₀ δ')
plt.ylabel('Desviación típica del remordimiento')
fig = plt.gcf()
fig.savefig('VarPseudoregUCB.pdf',format='pdf')
plt.show()