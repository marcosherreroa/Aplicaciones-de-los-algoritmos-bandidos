# -*- coding: utf-8 -*-
"""
Impacto Delta UCB Normal
Sección 8
Autor: Marcos Herrero
"""

import math
import random
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np

def samplePseudoRegretRandomChoice(n,k,arms,gaps):
    regret = 0
    
    for i in range(n):
        regret += gaps[random.choice(range(k))]
    
    return regret
    
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

n = 500
nDeltas = 40
samplenum = 200

regretRandom = np.zeros(nDeltas)
regretUCB = np.zeros(nDeltas)
C_UCB = np.zeros(nDeltas)
C_UCB2 = np.zeros(nDeltas)
Ctrivial = np.zeros(nDeltas)

arms = [None,None]
arms[1] = stats.norm(0,1)
gaps = [0,0]

Deltas = np.linspace(0.05,5,nDeltas)

for i in range(nDeltas):
    Delta = Deltas[i]
    
    arms[0] = stats.norm(Delta,1)
    gaps[1] = Delta
    
    for s in range(samplenum):
        regretRandom[i] += samplePseudoRegretRandomChoice(n, 2, arms, gaps)
        
    regretRandom[i] /= samplenum
    
    delta = 1/n**2
    for s in range(samplenum):
        regretUCB[i] += samplePseudoRegretUCB(n, 2, delta, arms, gaps)
    
    regretUCB[i]/= samplenum
    
    Ctrivial[i] = n*Delta  
    C_UCB[i] = 3*Delta+ 16*np.log(n)/Delta
    C_UCB2[i] = 8*np.sqrt(n*2*np.log(n))+3*Delta


fig = plt.figure()
plt.plot(Deltas,regretRandom,color= 'tab:green',label = 'Elección aleatoria')
plt.plot(Deltas,regretUCB, color = 'tab:blue', label='UCB (δ = 1/n²)')
plt.plot(Deltas,C_UCB,'--', color = 'tab:olive',label='C_UCB')
plt.plot(Deltas,C_UCB2,'--', color = 'tab:purple',label='C_UCB2')
plt.plot(Deltas,Ctrivial,'--',color = 'tab:red',label = 'n∆')
plt.xlabel('∆')
plt.ylabel('Remordimiento esperado')
plt.legend(loc = 'upper right')
fig.savefig('ImpactoDeltaNormalUCB1.pdf',format='pdf')
plt.show()

fig = plt.figure()
plt.plot(Deltas,regretRandom,color= 'tab:green',label = 'Elección aleatoria')
plt.plot(Deltas,regretUCB, color = 'tab:blue', label='UCB (δ = 1/n²)')
plt.plot(Deltas,C_UCB,'--', color = 'tab:olive',label='C_UCB')
#plt.plot(Deltas,C_UCB2,'--', color = 'tab:purple',label='C_UCB2')
#plt.plot(Deltas,Ctrivial,'--',color = 'tab:red',label = 'n∆')
plt.xlabel('∆')
plt.ylabel('Remordimiento esperado')
plt.ylim(0,100)
plt.legend(loc = 'upper right')
fig.savefig('ImpactoDeltaNormalUCB2.pdf',format='pdf')
plt.show()

arms = [None,None]
arms[1] = stats.norm(0,1)
gaps = [0,0]

invDeltas = np.arange(2,2+nDeltas)
regretRandom = np.zeros(nDeltas)
regretUCB = np.zeros(nDeltas)

for i in range(nDeltas):
    Delta = 1/(invDeltas[i])
    
    arms[0] = stats.norm(Delta,1)
    gaps[1] = Delta
    
    for s in range(samplenum):
        regretRandom[i] += samplePseudoRegretRandomChoice(n, 2, arms, gaps)
        
    regretRandom[i] /= samplenum
    
    delta = 1/n**2
    for s in range(samplenum):
        regretUCB[i] += samplePseudoRegretUCB(n, 2, delta, arms, gaps)
    
    regretUCB[i]/= samplenum
    
    Ctrivial[i] = n*Delta
    C_UCB[i] = 3*Delta+ 16*np.log(n)/Delta
    C_UCB2[i] = 8*np.sqrt(n*2*np.log(n))+3*Delta

fig = plt.figure()
plt.plot(invDeltas,regretRandom,color =  'tab:green',label='Elección aleatoria')
plt.plot(invDeltas,regretUCB, color = 'tab:blue', label='UCB (δ = 1/n²)')
plt.plot(invDeltas,C_UCB,'--',color = 'tab:olive',label ='C_UCB')
plt.plot(invDeltas,C_UCB2,'--',color ='tab:purple',label='C_UCB2')
plt.plot(invDeltas,Ctrivial,'--',color = 'tab:red',label ='n∆')
plt.xlabel('1/∆')
plt.ylabel('Remordimiento esperado')
plt.legend(loc='upper left')
fig.savefig('ImpactoDeltaNormalUCB3.pdf',format='pdf')
plt.show()

fig = plt.figure()
plt.plot(invDeltas,regretRandom,color =  'tab:green',label='Elección aleatoria')
plt.plot(invDeltas,regretUCB, color = 'tab:blue', label='UCB (δ = 1/n²)')
#plt.plot(invDeltas,C_UCB,'--',color = 'tab:olive',label ='C_UCB')
#.plot(invDeltas,C_UCB2,'--',color ='tab:purple',label='C_UCB2')
plt.plot(invDeltas,Ctrivial,'--',color = 'tab:red',label ='n∆')
plt.xlabel('1/∆')
plt.ylabel('Remordimiento esperado')
plt.ylim(0,40)
plt.legend(loc='upper left')
fig.savefig('ImpactoDeltaNormalUCB4.pdf',format='pdf')
plt.show()