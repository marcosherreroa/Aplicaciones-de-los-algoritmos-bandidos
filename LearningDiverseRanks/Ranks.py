# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 13:27:11 2022

@author: marco
"""

import math
import random
import abc

class User:
    
    def __init__(self, relevants):
        self.relevants = relevants
    
    def clicked(self,docs):
        for i in range(len(docs)):
            if docs[i] in self.relevants:
                return i
        
        return -1

class Bandit(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def selectArm(self):
        pass
    
    @abc.abstractmethod
    def reward(self,rwd):
        pass

class BanditEF(Bandit):
     
    def __init__(self,k,m):
        self.k = k
        self.chosenArm = -1
        self.m = m
        self.t = 0
        self.acum = k*[0]
        
    
    def selectArm(self):
        if self.t < self.m*self.k:
            self.chosenArm = self.t%self.k
        
        else:
            self.chosenArm = self.bestarm
        
        return self.chosenArm
    
    def reward(self,rwd):
        if self.t < self.m*self.k:
            self.acum[self.chosenArm] += rwd
        
        self.t += 1
        if self.t == self.m*self.k:
            print([i for i in range(self.k) if self.acum[i] == max(self.acum)])
            self.bestarm = random.choice([i for i in range(self.k) if self.acum[i] == max(self.acum)])


class BanditUCB(Bandit):
    
    def __init__(self,k,delta):
        self.k = k
        self.chosenArm = -1
        self.delta = delta
        self.t = 0
        self.meanReward = k*[0]
        self.T = k*[0]
        self.UCB = k*[math.inf]
        
    
    def selectArm(self):
        self.chosenArm = max(range(self.k),key=lambda i: self.UCB[i])
        return self.chosenArm
    
    def reward(self,rwd):
        self.meanReward[self.chosenArm] =  self.T[self.chosenArm]/(self.T[self.chosenArm]+1) \
                           *self.meanReward[self.chosenArm] + rwd/(self.T[self.chosenArm]+1)
    
        self.T[self.chosenArm] += 1
        self.UCB[self.chosenArm] = self.meanReward[self.chosenArm] + \
                                  math.sqrt((2*math.log(1/self.delta))/self.T[self.chosenArm])
            
        
        
numDocs = 100
p = 2
n = 10000

users = [User([1,3,4]),User([3,5,99]),User([5,60])]
indUsers = 0

bandits = [BanditUCB(numDocs, 0.3) for i in range(p)]


totalrwd = 0

for t in range(1,n):
    selected = p*[-1]
    presented = p*[-1]
    takeable = set(range(numDocs))
    
    for i in range(p):
        selected[i] = bandits[i].selectArm()
        if selected[i] not in takeable:
            presented[i] = random.choice(list(takeable))
        else:
            presented[i] = selected[i]
        takeable.remove(presented[i])
    
    print('Presentados: {}'.format(presented))
    
    indChosen = users[indUsers].clicked(presented)
    indUsers = (indUsers +1)%len(users)
    
    print('Elegido: {}'.format(indChosen))
    
    if indChosen != -1:
        totalrwd += 1
    
    for i in range(p):
        if i == indChosen and presented[i] == selected[i]:
            bandits[i].reward(1)
        
        else:
            bandits[i].reward(0)
    

print(totalrwd)
        