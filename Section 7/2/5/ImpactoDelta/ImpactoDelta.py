# -*- coding: utf-8 -*-

'''
Bandidos estocásticos: introducción, algoritmos y experimentos
TFG Informática
Sección 7.2.5
Figuras 8 y 9
Autor: Marcos Herrero Agustín
'''

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

n = 2000
samplenum = 500

arms = [None,None]
arms[1] = stats.uniform(0,1)
gaps = [0,0]

nDeltas = 100
invDeltas = np.arange(2,2+nDeltas)
regretRandom = np.zeros(nDeltas)
regretEF = np.zeros(nDeltas)
Ctrivial = np.empty(nDeltas) # n*Delta
C_EP2Min2arg = np.empty(nDeltas) # Delta + 4/Delta*max(...)

for i in range(nDeltas):
    Delta = 1/(invDeltas[i])
    
    arms[0] = stats.uniform(0,1+2*Delta)
    gaps[1] = Delta
    
    for s in range(samplenum):
        regretRandom[i] += samplePseudoRegretRandomChoice(n, 2, arms, gaps)
        
    regretRandom[i] /= samplenum
    
    m = computemTeor(n,Delta)
    for s in range(samplenum):
        regretEF[i] += samplePseudoRegretEF(n, 2, m, arms, gaps)
    
    regretEF[i]/= samplenum
    
    Ctrivial[i] = n*Delta
    
    if Delta == 0:
        C_EP2Min2arg[i] = 0
    else:
        C_EP2Min2arg[i] = Delta+4/Delta*(1+max(0,math.log(n*Delta*Delta/4)))
        

fig = plt.figure()
plt.plot(invDeltas,regretRandom,color =  'tab:green',label='Elección aleatoria')
plt.plot(invDeltas,regretEF, color = 'tab:blue', label='EP (m = m_Teor)')
plt.plot(invDeltas,C_EP2Min2arg,'--', color = 'tab:blue',label='C_2arg')
plt.plot(invDeltas,Ctrivial,'--',color = 'tab:red',label ='n∆')
plt.xlabel('1/∆')
plt.ylabel('Remordimiento esperado')
plt.legend(loc='upper left')
fig.savefig('ImpactoDelta3.pdf',format='pdf')
plt.show()

fig = plt.figure()
plt.plot(invDeltas[20:],regretRandom[20:],color =  'tab:green',label='Elección aleatoria')
plt.plot(invDeltas[20:],regretEF[20:], color = 'tab:blue', label='EP (m = m_Teor)')
#plt.plot(invDeltas[20:],C_EP2Min2arg[20:],'--', color = 'tab:blue',label='C_2arg')
plt.plot(invDeltas[20:],Ctrivial[20:],'--',color = 'tab:red',label ='n∆')
plt.xlabel('1/∆')
plt.ylabel('Remordimiento esperado')
#plt.ylim(0,20)
plt.legend(loc='upper left')
fig.savefig('ImpactoDelta4.pdf',format='pdf')
plt.show()
