# Bernoulli
# Aplicaciones de los algoritmos bandidos
# Autor : Marcos Herrero Agustín


import random
import sympy.stats as stats
import matplotlib.pyplot as plt


def randomPolicy(n,machines):
    chosen = n*[-1]
    totalrwd = 0
    
    for i in range(n):
        chosen[i] = random.randint(0,len(machines)-1)
        totalrwd += stats.sample(machines[chosen[i]])
    
    return chosen,totalrwd

def exploremtimes(n,machines,m):
    chosen = n*[-1]
    rwdbyArm = len(machines)*[0]
    timesArm = len(machines)*[0]
    totalrwd = 0
    ind = 0
    
    for i in range(m):
        for j in range(len(machines)):
            chosen[ind] = j
            rwd = stats.sample(machines[chosen[ind]])
            totalrwd += rwd
            rwdbyArm[j] += rwd
            timesArm[j] += 1
            ind += 1
    
    best = 0
    bestVal = rwdbyArm[0]/timesArm[0]
    for i in range(1,len(machines)):
        val = rwdbyArm[i]/timesArm[i]
        if val > bestVal:
            bestVal = val
            best = i
    
    print(best)    
    while ind < n:
        chosen[ind] = best
        rwd = stats.sample(machines[chosen[ind]])
        totalrwd += rwd
        ind += 1
    
    return chosen,totalrwd
    
def computeCumRegret(n,testnum,probs,chosen):
    
    regret = n*[0]
    acum = testnum*[0]
    maxprob = max(probs)
    for j in range(n):
        for i in range(testnum):
            #regret[j] += (maxprob - probs[chosen[i][j]])
            acum[i] += (maxprob - probs[chosen[i][j]])
            regret[j] += acum[i]
        
        regret[j] /= testnum
    
    return regret

#probs = [0.5 + i/20 for i in range(10)]   
#machines = 10*[None]
nmachines = 2
probs = [0.2,0.7]
machines = nmachines*[None]

for i in range(nmachines):
    
    pmf = {0: 1 - probs[i], 1: probs[i]}
    machines[i] = stats.FiniteRV('B('+str(i)+')',pmf)


n = 50
testnum = 100
chosen = testnum*[None]
totalrwd = testnum*[0]


for i in range(testnum):
    chosen[i], totalrwd[i] = randomPolicy(n,machines)


regret = computeCumRegret(n,testnum,probs,chosen)
plt.plot(regret,color='tab:blue')

for i in range(testnum):
    chosen[i], totalrwd[i] = exploremtimes(n,machines,1)

regret = computeCumRegret(n,testnum,probs,chosen)
plt.plot(regret,color='orange')


for i in range(testnum):
    chosen[i], totalrwd[i] = exploremtimes(n,machines,4)

regret = computeCumRegret(n,testnum,probs,chosen)
plt.plot(regret,color='green')

for i in range(testnum):
    chosen[i], totalrwd[i] = exploremtimes(n,machines,8)

regret = computeCumRegret(n,testnum,probs,chosen)
plt.plot(regret,color='red')

for i in range(testnum):
    chosen[i], totalrwd[i] = exploremtimes(n,machines,12)

regret = computeCumRegret(n,testnum,probs,chosen)
plt.plot(regret,color='tab:purple')


for i in range(testnum):
    chosen[i], totalrwd[i] = exploremtimes(n,machines,16)

regret = computeCumRegret(n,testnum,probs,chosen)
plt.plot(regret,color='brown')



plt.legend(['aleatorio','1 exploración','4 exploraciones','8 exploraciones','12 exploraciones','16 exploraciones'])
plt.xlabel("Ronda")
plt.ylabel("Regret acumulado")
fig = plt.gcf()
fig.savefig("Bernoulli.pdf",format='pdf')
plt.show()
    

