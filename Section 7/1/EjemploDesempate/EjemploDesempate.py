# -*- coding: utf-8 -*-

'''
Bandidos estocásticos: introducción, algoritmos y experimentos
TFG Informática
Sección 7.1
Figuras 3, 4 y 5
Autor: Marcos Herrero Agustín
'''

import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np

def regretEFFavsFirst(n,m,arms,gaps):
    rwds = 2*[0]
    
    for i in range(m):
        rwds[0] += arms[0].rvs()
        rwds[1] += arms[1].rvs()
    
    bestarm = min([i for i in range(2) if rwds[i] == max(rwds)])
    
    return m*sum(gaps)+(n-2*m)*gaps[bestarm]

def regretEFFavsSecond(n,m,arms,gaps):
    rwds = 2*[0]
    
    for i in range(m):
        rwds[0] += arms[0].rvs()
        rwds[1] += arms[1].rvs()
    
    bestarm = max([i for i in range(2) if rwds[i] == max(rwds)])
    
    return m*sum(gaps)+(n-2*m)*gaps[bestarm]

n = 1000
nsamples = 100

ms = np.array(range(0,201,2))
numms = len(ms)

arms = 2*[None]
arms[0] = stats.norm(0,1)
arms[1] = stats.norm(-0.1,1)

gaps = 2*[0]
gaps[0] = 0
gaps[1] = 0.1


regretFavsFirst = numms*[0]
regretFavsSecond = numms*[0]

for i in range(numms):
    
    for k in range(nsamples):
        
        regretFavsFirst[i] = k/(k+1)*regretFavsFirst[i] + 1/(k+1)*regretEFFavsFirst(n,ms[i],arms,gaps)
        regretFavsSecond[i] = k/(k+1)*regretFavsSecond[i] + 1/(k+1)*regretEFFavsSecond(n,ms[i],arms,gaps)

colors = ['tab:blue','tab:red']
legend = ['Fav. brazo 1','Fav. brazo 2']

'''



box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(legend, loc = 'center left', bbox_to_anchor=(1, 0.5))                  
                 
fig.savefig('EFDesempate.pdf',format='pdf')
plt.show()
'''
 
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(ms,regretFavsFirst,color = colors[0])
ax.plot(ms,regretFavsSecond, color = colors[1])       
plt.xlabel('m')
plt.ylabel('Remordimiento esperado')
plt.legend(legend,loc = 'upper right')
fig.savefig('EFDesempate.pdf',format='pdf')
plt.show()

