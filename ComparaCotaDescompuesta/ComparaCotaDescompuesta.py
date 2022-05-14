# -*- coding: utf-8 -*-
"""
Compara Cota Descompuesta
Explora-Primero
Autor: Marcos Herrero
"""

import matplotlib.pyplot as plt
import numpy as np


C = 4*np.exp(-0.5)

#Para n fijo
n = 1000
nDeltas = 40

Ctrivial = np.empty(nDeltas) # n*Delta
C_EP2Min2arg = np.empty(nDeltas) # Delta + 4/Delta*max(...)
CDesc = np.empty(nDeltas) # Delta + C√n

Deltas = np.linspace(0,5,nDeltas)

for i in range(nDeltas):
    Delta = Deltas[i]
    
    Ctrivial[i] = n*Delta
    CDesc[i] = Delta + C*np.sqrt(n)
    
    if Delta == 0:
        C_EP2Min2arg[i] = 0
    else:
        C_EP2Min2arg[i] = Delta+4/Delta*(1+max(0,np.log(n*Delta*Delta/4)))

fig = plt.figure()
plt.plot(Deltas,C_EP2Min2arg,'--', color = 'tab:blue',label='C_2arg')
plt.plot(Deltas,Ctrivial,'--',color = 'tab:red',label ='n∆')
plt.plot(Deltas,CDesc,'--',color = 'tab:green',label = 'C_Desc')
plt.xlabel('∆')
plt.ylabel('Remordimiento esperado')
plt.ylim(0,100)
plt.legend(loc = 'upper right')
fig.savefig('CompDesc1.pdf',format='pdf')
plt.show()


invDeltas = np.arange(1,nDeltas+1)

for i in range(nDeltas):
    Delta = 1/invDeltas[i]
    
    Ctrivial[i] = n*Delta
    CDesc[i] = Delta + C*np.sqrt(n)
    
    if Delta == 0:
        C_EP2Min2arg[i] = 0
    else:
        C_EP2Min2arg[i] = Delta+4/Delta*(1+max(0,np.log(n*Delta*Delta/4)))

fig = plt.figure()
plt.plot(invDeltas,C_EP2Min2arg,'--', color = 'tab:blue',label='C_2arg')
plt.plot(invDeltas,Ctrivial,'--',color = 'tab:red',label ='n∆')
plt.plot(invDeltas,CDesc,'--',color = 'tab:green',label = 'C_Desc')
plt.xlabel('1/∆')
plt.ylabel('Remordimiento esperado')
plt.ylim(0,100)
plt.legend(loc = 'upper right')
fig.savefig('CompDesc2.pdf',format='pdf')
plt.show()

# Para Delta fijo
Delta = 1
nns = 40

Ctrivial = np.empty(nns) # n*Delta
C_EP2Min2arg = np.empty(nns) # Delta + 4/Delta*max(...)
CDesc = np.empty(nns) # Delta + C√n

ns = np.arange(1,10*nns,10)

for i in range(nns):
    n = ns[i]
    
    Ctrivial[i] = n*Delta
    CDesc[i] = Delta + C*np.sqrt(n)
    
    if Delta == 0:
        C_EP2Min2arg[i] = 0
    else:
        C_EP2Min2arg[i] = Delta+4/Delta*(1+max(0,np.log(n*Delta*Delta/4)))

fig = plt.figure()
plt.plot(ns,C_EP2Min2arg,'--', color = 'tab:blue',label='C_2arg')
plt.plot(ns,Ctrivial,'--',color = 'tab:red',label ='n∆')
plt.plot(ns,CDesc,'--',color = 'tab:green',label = 'C_Desc')
plt.xlabel('n')
plt.ylabel('Remordimiento esperado')
plt.ylim(0,100)
plt.legend(loc = 'upper right')
fig.savefig('CompDesc3.pdf',format='pdf')
plt.show()

Delta = 0.01
nns = 40

Ctrivial = np.empty(nns) # n*Delta
C_EP2Min2arg = np.empty(nns) # Delta + 4/Delta*max(...)
CDesc = np.empty(nns) # Delta + C√n

ns = np.arange(1,10*nns,10)

for i in range(nns):
    n = ns[i]
    
    Ctrivial[i] = n*Delta
    CDesc[i] = Delta + C*np.sqrt(n)
    
    if Delta == 0:
        C_EP2Min2arg[i] = 0
    else:
        C_EP2Min2arg[i] = Delta+4/Delta*(1+max(0,np.log(n*Delta*Delta/4)))

fig = plt.figure()
plt.plot(ns,C_EP2Min2arg,'--', color = 'tab:blue',label='C_2arg')
plt.plot(ns,Ctrivial,'--',color = 'tab:red',label ='n∆')
plt.plot(ns,CDesc,'--',color = 'tab:green',label = 'C_Desc')
plt.xlabel('n')
plt.ylabel('Remordimiento esperado')
plt.ylim(0,100)
plt.legend(loc = 'upper right')
fig.savefig('CompDesc4.pdf',format='pdf')
plt.show()

