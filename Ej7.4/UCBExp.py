# -*- coding: utf-8 -*-
"""
Ejercicio 7.8 de Tor Lattimore
Autor: Marcos Herrero
"""

import math
import random
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np

def computemTeor(n,Delta):
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
    
        
def plotDeltaRegret():
       
    n = 1000
    sampleNum = 300
    
    arms = 2*[None]
    arms[0] =  stats.norm(0,1)
    gaps = 2*[0]
    
    nDeltas = 20
    Deltas = np.linspace(0,1,nDeltas)
    regretEF25 = np.empty(nDeltas)
    regretEF50 = np.empty(nDeltas)
    regretEF75 = np.empty(nDeltas)
    regretEF100 = np.empty(nDeltas)
    regretEFmTeor = np.empty(nDeltas)
    regretEFOptimo = np.empty(nDeltas)
    regretUCB = np.empty(nDeltas)
    
    mTeor = nDeltas*[0]
    mOpt = nDeltas*[0]
    
    for i in range(nDeltas):
        Delta = Deltas[i]
        
        arms[1]= stats.norm(-Delta,1)
        gaps[1] = Delta
        
        regretEF25[i] = 0
        for k in range(sampleNum):
            regretEF25[i] += samplePseudoRegretEF(n,2,25,arms,gaps)
        
        regretEF25[i] /= sampleNum
        
        regretEF50[i] = 0
        for k in range(sampleNum):
            regretEF50[i] += samplePseudoRegretEF(n,2,50,arms,gaps)
        
        regretEF50[i] /= sampleNum
        
        regretEF75[i] = 0
        for k in range(sampleNum):
            regretEF75[i] += samplePseudoRegretEF(n,2,75,arms,gaps)
        
        regretEF75[i] /= sampleNum
        
        regretEF100[i] = 0
        for k in range(sampleNum):
            regretEF100[i] += samplePseudoRegretEF(n,2,100,arms,gaps)
        
        regretEF100[i] /= sampleNum
        
        regretEFmTeor[i]= 0
        mTeor[i] = computemTeor(n,Delta)
        for k in range(sampleNum):
            regretEFmTeor[i] += samplePseudoRegretEF(n,2,mTeor[i],arms,gaps)
        
        regretEFmTeor[i] /= sampleNum
        
        regretEFOptimo[i] = 0
        mOpt[i] = computemOpt(n,Delta) 
        for k in range(sampleNum):
            regretEFOptimo[i] += samplePseudoRegretEF(n,2,mOpt[i],arms,gaps)
        
        regretEFOptimo[i] /= sampleNum
        
        regretUCB[i] = 0
        for k in range(sampleNum):
            regretUCB[i] += samplePseudoRegretUCB(n,2,1/(n*n),arms,gaps)
        
        regretUCB[i] /= sampleNum
    
    fig = plt.figure()
    plt.plot(Deltas,regretEF25, color='tab:blue',label= 'EF (m=25)')
    plt.plot(Deltas,regretEF50, color='tab:green',label = 'EF (m=50)')
    plt.plot(Deltas,regretEF75, color='tab:olive',label = 'EF (m=75)')
    plt.plot(Deltas,regretEF100, color='tab:red', label = 'EF (m=100)')
    plt.plot(Deltas,regretEFmTeor, color='tab:purple',label = 'EF (m=m_Teor)')
    plt.plot(Deltas,regretEFOptimo, color='tab:gray', label = 'EF (m=m_Opt)')
    plt.plot(Deltas,regretUCB, color='black', label = 'UCB')
    plt.xlabel('∆')
    plt.ylabel('Remordimiento esperado')
    plt.legend(loc='upper left',ncol = 2)
    fig.savefig('UCBDeltaRegret.pdf',format='pdf')
    plt.show()
    
    fig = plt.figure()
    plt.plot(Deltas, mTeor, color='tab:purple', label = 'm_Teor')
    plt.plot(Deltas,mOpt, color = 'tab:gray', label = 'm_Opt')
    plt.xlabel('∆')
    plt.ylabel('m')
    plt.legend(loc='upper left')
    fig.savefig('ms.pdf',format='pdf')
    plt.show()
    

plotDeltaRegret()

