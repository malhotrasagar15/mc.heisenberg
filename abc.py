# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 10:54:30 2021

@author: Aman Arya
"""
import random
import math
import numpy as np
from matplotlib import pyplot as plt

def func():
        theta = np.random.random()*np.pi
        phi = np.random.random()*2*np.pi
        p = math.sin(theta) * math.cos(phi)
        q = math.sin(theta) * math.sin(phi)
        r = math.cos(theta)
        return p,q,r
func()
n=10
lattice = [[[ func() for i in range(n) ] for j in range(n) ] for k in range(n)]

def ising_model(n,tem,simulaions,J):      #tem=temperature
    energy_tem = []                               
    n_cube = n*n*n
    
    for t in tem:                  
        energy = []
        energy_sq = []
        
        size = len(lattice[0])
       
        for simulation in range(simulaions): 
            for l in range(n_cube):           #Each simulation is (n^2) spin flips acc to mcmc rule
                x = random.randint(0,n-1)
                y = random.randint(0,n-1)
                z = random.randint(0,n-1)
                                                            #cal. change in energy
                DeltaE = 2*lattice[x,y,z]*J*(lattice[x,(y-1)%size,z] + lattice[x,(y+1)%size,z] + lattice[(x-1)%size,y,z] + lattice[(x+1)%size,y,z] + lattice[x,y,(z+1)%size] + lattice[x,y,(z-1)%size])
                if DeltaE<=0:                                    # here metropolis alg.
                    lattice[x,y,z]=lattice[x,y,z]*-1       
                elif random.random()<(np.exp(-DeltaE/t)):   
                    lattice[x,y,z]=lattice[x,y,z]*-1 
           
            energyy = 0   
            for i in range(size):                    #cal variables after spin flip
                for j in range(size):
                    for k in range(size):
                        energyy=energyy+ lattice[i,j,k]*(lattice[i,(j-1)%size,k] + lattice[i,(j+1)%size,k] + lattice[(i-1)%size,j,k] + lattice[(i+1)%size,j,k] + lattice[i,j,(k+1)%size] + lattice[i,j,(k-1)%size])
            energy.append((-J*energyy)/n_cube)
            energy_sq.append(((-J*energyy)/n_cube)**2)  
            
        
        energy_tem.append(np.mean(energy[200:]))   #append the average energy per site for all simulations
    return (tem,energy_tem)
tem = np.linspace(0.50,10.0,10)
X = ising_model(16,tem,1000,1)

plt.plot(X[0],X[1])
plt.xlabel('temerature (K)')
plt.ylabel('Energy per site')
