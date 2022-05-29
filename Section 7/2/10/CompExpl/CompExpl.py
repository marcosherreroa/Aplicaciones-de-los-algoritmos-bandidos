# -*- coding: utf-8 -*-

"""
Bandidos estocásticos: introducción, algoritmos y experimentos
TFG Informática
Sección 7.2.10
Figuras 15 y 16
Autor: Marcos Herrero Agustín
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
    
def samplePseudoRegretRandomChoice(n,k,arms,gaps):
    regret = 0
    
    for i in range(n):
        regret += gaps[random.choice(range(k))]
    
    return regret

def samplePseudoRegretEF(n,k,m,arms,gaps):
    rwds = k*[0]
    
    for i in range(m):
        for j in range(k):
            rwds[j] += arms[j].rvs()
    
    maximum = max(rwds)
    bestarm = random.choice([i for i in range(k) if rwds[i] == maximum])
    
    return m*sum(gaps)+(n-m*k)*gaps[bestarm]

n = 1000
sampleNum = 1000
 
arms = 2*[None]
arms[0] =  stats.bernoulli(0.95)
gaps = 2*[0]

nDeltas = 40
regretmExpl = np.zeros(nDeltas)
regretmTeor = np.zeros(nDeltas)
regretRandom = np.zeros(nDeltas)
C_EP2Min = np.empty(nDeltas)

mExpl = computemExpl(n)

mExpl = nDeltas*[mExpl]  
mTeor = nDeltas*[0]
mOpt = nDeltas*[0]

#Para delta grande
Deltas = np.linspace(0.05,0.9,nDeltas)

for i in range(nDeltas):
    Delta = Deltas[i]
    
    arms[1]= stats.bernoulli(0.95 - Delta)
    gaps[1] = Delta
    
    mTeor[i] = computemTeor(n,Delta)
    
    for k in range(sampleNum):
        regretmExpl[i] += samplePseudoRegretEF(n,2,mExpl[i],arms,gaps)
        regretmTeor[i] += samplePseudoRegretEF(n,2,mTeor[i],arms,gaps)
        regretRandom[i] += samplePseudoRegretRandomChoice(n, 2, arms, gaps)
        
    
    regretmExpl[i] /= sampleNum
    regretmTeor[i] /= sampleNum
    regretRandom[i] /= sampleNum
    
    if Delta == 0:
        C_EP2Min[i] = 0
    else:
        C_EP2Min[i] = min(n*Delta,Delta+4/Delta*(1+max(0,math.log(n*Delta*Delta/4))))

 

fig = plt.figure()
plt.plot(Deltas,regretRandom, color='tab:green', label = 'Elección aleatoria')
plt.plot(Deltas,regretmExpl,color='tab:blue',label='EP (m=m_Expl)')
plt.plot(Deltas,regretmTeor, color='tab:purple',label = 'EP (m = m_Teor)')
plt.plot(Deltas,C_EP2Min,'--', color='tab:purple',label= 'C_EP2_min')
plt.xlabel('∆')
plt.ylabel('Remordimiento esperado')
plt.ylim(0,80)
plt.legend(loc='upper right')
fig.savefig('RegretExpl1.pdf',format='pdf')
plt.show()

fig = plt.figure()
plt.plot(Deltas, mExpl, color='tab:blue', label = 'm_Expl')
plt.plot(Deltas, mTeor, color='tab:purple', label = 'm_Teor')
plt.xlabel('∆')
plt.ylabel('m')
plt.legend(loc='upper right')
fig.savefig('msExpl1.pdf',format='pdf')
plt.show()

#Para Delta pequeño
invDeltas = np.arange(2,2+nDeltas)

for i in range(nDeltas):
    Delta = 1/invDeltas[i]
    
    arms[1]= stats.bernoulli(0.95 - Delta)
    gaps[1] = Delta
    
    mTeor[i] = computemTeor(n,Delta)
    
    for k in range(sampleNum):
        regretmExpl[i] += samplePseudoRegretEF(n,2,mExpl[i],arms,gaps)
        regretmTeor[i] += samplePseudoRegretEF(n,2,mTeor[i],arms,gaps)
        regretRandom[i] += samplePseudoRegretRandomChoice(n, 2, arms, gaps)
        
    
    regretmExpl[i] /= sampleNum
    regretmTeor[i] /= sampleNum
    regretRandom[i] /= sampleNum
    
    if Delta == 0:
        C_EP2Min[i] = 0
    else:
        C_EP2Min[i] = min(n*Delta,Delta+4/Delta*(1+max(0,math.log(n*Delta*Delta/4))))

 

fig = plt.figure()
plt.plot(invDeltas,regretRandom, color='tab:green', label = 'Elección aleatoria')
plt.plot(invDeltas,regretmExpl,color='tab:blue',label='EP (m=m_Expl)')
plt.plot(invDeltas,regretmTeor, color='tab:purple',label = 'EP (m = m_Teor)')
plt.plot(invDeltas,C_EP2Min,'--', color='tab:purple',label= 'C_EP2_min')
plt.xlabel('1/∆')
plt.ylabel('Remordimiento esperado')
plt.ylim(0,80)
plt.legend(loc='upper right')
fig.savefig('RegretExpl2.pdf',format='pdf')
plt.show()

fig = plt.figure()
plt.plot(invDeltas, mExpl, color='tab:blue', label = 'm_Expl')
plt.plot(invDeltas, mTeor, color='tab:purple', label = 'm_Teor')
plt.xlabel('1/∆')
plt.ylabel('m')
plt.legend(loc='upper right')
fig.savefig('msExpl2.pdf',format='pdf')
plt.show()