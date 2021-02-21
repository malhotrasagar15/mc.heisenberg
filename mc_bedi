
%matplotlib inline
import random
import numpy as np
from matplotlib import pyplot as plt
import time

def lattice_initial(n):                 #function which gives the starting lattice state of size n X n
    return(2*np.random.randint(2, size=(n,n))-1)


def ising_model(n,tem,simulaions,J):  
    energy_tem = []         
    heat_cap_tem = [] 
    mag_tem = []
    susc_tem = []
    n_sq = n*n
    for t in tem:          
        start_tem = time.time()
        energy = []
        mag = []
        energy_sq = []
        mag_sq = []
        lattice = lattice_initial(n)   #initializing the lattice
        size = len(lattice[0])
        for simulation in range(simulaions): 
            for k in range(n_sq):           #Each simulation is (n^2) spin flips acc to mcmc rule
                x = random.randint(0,n-1)
                y = random.randint(0,n-1) 
                DeltaE = 2*lattice[x,y]*J*(lattice[x,(y-1)%size] + lattice[x,(y+1)%size] + lattice[(x-1)%size,y] + lattice[(x+1)%size,y])
                if DeltaE<=0: #when the flip lowers the energy of lattice, it is laways accepted. 
                    lattice[x,y]=lattice[x,y]*-1       
                elif random.random()<(np.exp(-DeltaE/t)):     #when the flip increases the energy of lattice, it is accepted with probability exp(-DeltaE/T).
                    lattice[x,y]=lattice[x,y]*-1 
            energyy = 0   
            for i in range(size):
                for j in range(size):
                    energyy=energyy+ lattice[i,j]*(lattice[i,(j-1)%size] + lattice[i,(j+1)%size] + lattice[(i-1)%size,j] + lattice[(i+1)%size,j])
            energy.append(1.0*(-J*energyy/2.0)/n_sq)
            energy_sq.append((1.0*(-J*energyy/2.0)**2)/n_sq**2)  
            mag.append(abs(sum(sum(lattice)))/(1.0*n_sq))
            mag_sq.append((abs(sum(sum(lattice)))/(1.0*n_sq))**2)
            
        plt.plot(range(simulaions),energy)
        plt.title("Energy of the lattice vs no. of iterations at tem "+str(t)+"K")
        plt.xlabel("No. of simulaions with each simulation having "+str(n*n)+" spin flips)")
        plt.ylabel("Energy per unit site")
        plt.show()
        energy_tem.append(np.mean(energy[200:])) #append the average energy per site for all simulations
        cv = (1.0/(t**2))*(np.mean(energy_sq[200:])-np.mean(energy[100:])**2) 
        heat_cap_tem.append(cv)
        mag_tem.append(np.mean(mag))  #append the average magnetization of lattice
        susc = (1.0/t)*(np.mean(mag_sq[200:])-np.mean(mag[200:])**2)
        susc_tem.append(susc)
        end_tem = time.time()
    return (tem,energy_tem,heat_cap_tem,mag_tem,susc_tem) 
tem = np.linspace(1.9,4.0,20)
X = ising_model(16,tem,1000,1)
plt.plot(X[0],X[2])
plt.xlabel('temerature (K)')
plt.ylabel('Heat Capacity (Cv)')
plt.show()
plt.plot(X[0],X[1])
plt.xlabel('temerature (K)')
plt.ylabel('Energy per site')
plt.show()
plt.plot(X[0],X[3])
plt.xlabel('temerature (K)')
plt.ylabel('Magnetization per spin')
plt.show()
plt.plot(X[0],X[4])
plt.xlabel('temerature (K)')
plt.ylabel('Susceptibility (Chi)')
