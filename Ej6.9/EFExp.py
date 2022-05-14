# -*- coding: utf-8 -*-
"""
Experimento ETC
Autor : Marcos Herrero
"""

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

def plotDeltaRegret():
       
    n = 1000
    sampleNum = 1000
    
    arms = 2*[None]
    arms[0] =  stats.norm(0,1)
    gaps = 2*[0]
    
    nDeltas = 40
    Deltas = np.linspace(0,1,nDeltas)
    regretmExpl = np.zeros(nDeltas)
    regretmTeor = np.zeros(nDeltas)
    regretmOpt = np.zeros(nDeltas)
    C_EP2Min = np.empty(nDeltas)
    C_EP3Min = np.empty(nDeltas)
    
    mExpl = computemExpl(n)
    
    mExpl = nDeltas*[mExpl]  
    mTeor = nDeltas*[0]
    mOpt = nDeltas*[0]
    
    for i in range(nDeltas):
        Delta = Deltas[i]
        
        arms[1]= stats.norm(-Delta,1)
        gaps[1] = Delta
        
        mTeor[i] = computemTeor(n,Delta)
        mOpt[i] = computemOpt(n,Delta)
        
        for k in range(sampleNum):
            regretmExpl[i] += samplePseudoRegretEF(n,2,mExpl[i],arms,gaps)
            regretmTeor[i] += samplePseudoRegretEF(n,2,mTeor[i],arms,gaps)
            regretmOpt[i] += samplePseudoRegretEF(n,2,mOpt[i],arms,gaps)
            
        
        regretmExpl[i] /= sampleNum
        regretmTeor[i] /= sampleNum
        regretmOpt[i] /= sampleNum
        
        if Delta == 0:
            C_EP2Min[i] = 0
        else:
            C_EP2Min[i] = min(n*Delta,Delta+4/Delta*(1+max(0,math.log(n*Delta*Delta/4))))
        
        C_EP3Min[i] = (Delta+math.sqrt(2)*math.exp(-1/2))*math.ceil(n**(2/3))
    
    fig = plt.figure()
    plt.plot(Deltas,regretmExpl,color='tab:red',label='EP (m=m_Expl)')
    plt.plot(Deltas,regretmTeor, color='tab:green',label = 'EP (m = m_Teor)')
    plt.plot(Deltas,regretmOpt, color='tab:blue', label = 'EP (m = m_Opt)')
    plt.plot(Deltas,C_EP2Min,'--', color='tab:green',label= 'C_EP2_min')
    plt.plot(Deltas,C_EP3Min,'--',color='tab:red',label='C_EP3_min')
    plt.xlabel('∆')
    plt.ylabel('Remordimiento esperado')
    plt.ylim(0,80)
    plt.legend(loc='upper right')
    fig.savefig('EFExpDeltaRegret.pdf',format='pdf')
    plt.show()
    
    fig = plt.figure()
    plt.plot(Deltas, mExpl, color='tab:red', label = 'm_Expl')
    plt.plot(Deltas, mTeor, color='tab:green', label = 'm_Teor')
    plt.plot(Deltas,mOpt, color = 'tab:blue', label = 'm_Opt')
    plt.xlabel('∆')
    plt.ylabel('m')
    plt.legend(loc='upper right')
    fig.savefig('ms.pdf',format='pdf')
    plt.show()

def plotmRegretAndVar():
    
    n = 1000
    Delta = 0.1
    sampleNum = 600
    
    arms = 2*[None]
    arms[0] = stats.Normal('X0',0,1)
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
    
plotDeltaRegret()
#plotmRegretAndVar()