import random
import math
import numpy as np
from operator import add
from matplotlib import pyplot as plt

def func():
        theta = np.random.random()*np.pi
        phi = np.random.random()*2*np.pi
        p = math.sin(theta) * math.cos(phi)
        q = math.sin(theta) * math.sin(phi)
        r = math.cos(theta)
        return [p,q,r]

n=10
'''def lattice_initial(n):                 #function which gives the starting lattice state of size n X n
    return(2*np.random.randint(2, size=(n,n,n))-1)'''


lattice = [[[ func() for i in range(n) ] for j in range(n) ] for k in range(n)]

'''Z = [lattice[0][1][0],lattice[1][0][1]]
print(len(Z))'''

def addmul(L):
        l=np.zeros(3)
        print('l',np.shape(l))
        print('L',np.shape(L[0]))
        
        
        for i in range(0,len(L)):
                l = np.add(l,L[i])
        return l

'''Z = [lattice[0][1][0],lattice[1][0][1],lattice[0][0][0]]
print(Z)
print(np.add(Z[0],Z[1]))
print(addmul(Z))
print(np.dot(Z[0],addmul(Z)))'''


'''print(lattice)
print(lattice[0][0][0],lattice[0][1][0],lattice[0][1][1])
print(np.add(lattice[0][0][0],np.add(lattice[0][1][0],lattice[0][1][1])))'''

def ising_model(n,tem,simulaions,J):      
    energy_tem = []                               
    n_cube = n*n*n
    
    for t in tem:                  
        energy = []
        energy_sq = []
        
           
        size = n
        
       
        for simulation in range(simulaions): 
            for c in range(n_cube):           
                x = random.randint(0,n-1)
                y = random.randint(0,n-1)
                z = random.randint(0,n-1)
                N = [lattice[x][(y-1)%size][z],lattice[x][(y+1)%size][z],lattice[(x-1)%size][y][z],lattice[(x+1)%size][y][z],lattice[x][y][(z+1)%size],lattice[x][y][(z-1)%size]]                                            
                M = [lattice[0][1][2],lattice[0][2][1],lattice[0][2][3],lattice[4][1][3],lattice[0][5][2],lattice[1][5][3]]
                DeltaE = 2*J*np.dot(lattice[x][y][z],addmul(N))
                if DeltaE<=0:                                    
                    lattice[x][y][z]=lattice[x][y][z]*-1       
                elif random.random()<(np.exp(-DeltaE/t)):   
                    lattice[x][y][z]=lattice[x][y][z]*-1 
           
            energyy = 0   
            for i in range(size):                    
                for j in range(size):
                    for k in range(size):
                        energyy=energyy+ np.dot(lattice[i][j][k],addmul([lattice[i][(j-1)%size][k],lattice[i][(j+1)%size][k],lattice[(i-1)%size][j][k],lattice[(i+1)%size][j][k],lattice[i][j][(k+1)%size],lattice[i][j][(k-1)%size]]))
            energy.append((-J*energyy)/n_cube)
            energy_sq.append(((-J*energyy)/n_cube)**2)
          
        
        energy_tem.append(np.mean(energy[200:]))   
    return (tem,energy_tem)
tem = np.linspace(0.50,10.0,10)
X = ising_model(10,tem,1000,1)

plt.plot(X[0],X[1])
plt.xlabel('temerature (K)')
plt.ylabel('Energy per site')
