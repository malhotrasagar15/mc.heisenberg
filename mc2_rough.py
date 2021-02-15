%matplotlib inline
import random
import math
import numpy as np
from matplotlib import pyplot as plt
import time

def latt():
    b = math.pi
    def theta(b):
        return b*random.randint(0,10.0)/10.0
    def phi(b):
        return b*random.randint(0,10.0)/5.0
    x=(math.sin(theta(b)))*(math.cos(phi(b)))
    y=(math.sin(theta(b)))*(math.sin(phi(b)))
    z=(math.cos(theta(b)))

n=10
a = np.zeros((n,n,n,2))
for i in range(n):
    for j in range(n):
        for k in range(n):
            a[i][j][k] = latt()


def ising_model(n,tem,simulaions,J):      #tem=temperature
    energy_tem = []                               
    ncube = n*n*n
    
    for t in tem:          
        start_tem = time.time()      #start of simulation
        
        energy = []
        energy_sq = []
        
        lattice = latt(n)   #initializing the lattice
        size = len(lattice[0])
       
        for simulation in range(simulaions): 
            for kk in range(ncube):           #Each simulation is (n^2) spin flips acc to mcmc rule
                x = random.randint(0,n-1)
                y = random.randint(0,n-1)
                z = random.randint(0,n-1)
                                                            #cal. change in energy
                DeltaE = 2*lattice[x,y]*J*(lattice[x,(y-1)%size,z] + lattice[x,(y+1)%size,z] + lattice[(x-1)%size,y,z] + lattice[(x+1)%size,y,z] + lattice[x,y,(z+1)%size] + lattice[x,y,(z-1)%size])
                if DeltaE<=0:                                    # here metropolis alg.
                    lattice[x,y,z]=lattice[x,y,z]*-1       
                elif random.random()<(np.exp(-DeltaE/t)):   
                    lattice[x,y,z]=lattice[x,y]*-1 
            energyy = 0   
            for i in range(size):                    #cal variables after spin flip
                for j in range(size):
                    for k in range(size):
                        
                        energyy=energyy+ lattice[i,j]*(lattice[i,(j-1)%size] + lattice[i,(j+1)%size] + lattice[(i-1)%size,j] + lattice[(i+1)%size,j])
            energy.append((-J*energyy)/ncube)
            energy_sq.append(((-J*energyy)/ncube)**2)  
            
        
        energy_tem.append(np.mean(energy[200:]))   #append the average energy per site for all simulations
        end_tem = time.time()
    return (tem,energy_tem)
tem = np.linspace(0.50,10.0,21)
X = ising_model(16,tem,1000,1)

plt.plot(X[0],X[1])
plt.xlabel('temerature (K)')
plt.ylabel('Energy per site')
