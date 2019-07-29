#!/usr/bin/env python
"""
Calculate 2D Ising model using Metropolis Monte-Carlo method.
*site fixed*
Depends on numpy, metaplotlib, tqdm and mpi4py
"""

#Import library
from __future__ import print_function
from __future__ import division
import numpy as np
from numpy.random import rand
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Initialize MPI part
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

#----------------------------------------------------------------------
##  BLOCK OF FUNCTIONS USED IN THE MAIN CODE
#----------------------------------------------------------------------
def read_input():
    ''' a subroutine to get information from input.MC'''
    dataset={}
    file = open('input.MC', "r")
    #Default
    data = file.readlines()
    nt      = 18         #  number of temperature points
    N       = 16         #  size of the lattice, N x N
    eqSteps = 8000       #  number of MC sweeps for equilibration
    mcSteps = 4000       #  number of MC sweeps for calculation
    D_data =  1.
    T_low      =  1.53
    T_high     =  3.28

    for line in data:
        key, value = line.split("=")
        dataset[key.strip()] = value.strip()
    # Read data
    nt          = int(dataset["nt"])
    N           = int(dataset["N"])
    eqSteps     = int(dataset["eqSteps"])
    mcSteps     = int(dataset["mcSteps"])
    D_data      = float(dataset["D_data"])
    T_low       = float(dataset["T_low"])
    T_high      = float(dataset["T_high"])
    return nt, N, eqSteps, mcSteps, D_data, T_low, T_high

def initialstate(N):
    ''' generates a random spin configuration for initial condition'''
    state = 2*np.random.randint(2, size=(N,N))-1
    return state


def mcmove(config, beta,D_data):
    '''Monte Carlo move using Metropolis algorithm '''
    for i in range(N):
        for j in range(N):
                a = np.random.randint(0, N)
                b = np.random.randint(0, N)
                s =  config[a, b]
                nb = config[(a+1)%N,b] + config[a,(b+1)%N] + config[(a-1)%N,b] + config[a,(b-1)%N]
                nnb = config[(a+1)%N,(b+1)%N] + config[(a+1)%N,(b-1)%N] + config[(a-1)%N,(b-1)%N] + config[(a-1)%N,(b+1)%N]
                cost = 2*s*D_data*nb#+2*0.5*s*nnb
                if cost < 0:
                    s *= -1
                elif rand() < np.exp(-cost*beta):
                    s *= -1
                config[a, b] = s
    return config


def calcEnergy(config):
    '''Energy of a given configuration'''
    energy = 0
    for i in range(len(config)):
        for j in range(len(config)):
            S = config[i,j]
            nb = config[(i+1)%N, j] + config[i,(j+1)%N] + config[(i-1)%N, j] + config[i,(j-1)%N]
            energy += -nb*S
    return energy/4.


def calcMag(config):
    '''Magnetization of a given configuration'''
    mag = np.abs(np.sum(config))
    return mag

#----------------------------------------------------------------------
#  Initialize Model parameters
#----------------------------------------------------------------------
# change these parameters for a smaller (faster) simulation
nt, N, eqSteps, mcSteps, D_data, T_low, T_high = read_input()
#----------------------------------------------------------------------
T       = np.linspace(T_low, T_high, nt);
E,M,C,X = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
E_1,M_1,C_1,X_1 = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
n1, n2  = 1.0/(mcSteps*N*N), 1.0/(mcSteps*mcSteps*N*N)
# divide by number of samples, and by system size to get intensive values

#----------------------------------------------------------------------
#  MAIN PART OF THE CODE
#----------------------------------------------------------------------

for tt in range(int(nt*rank/size),int(nt*(rank+1)/size)):
    E1 = M1 = E2 = M2 = 0
    config = initialstate(N)
    iT=1.0/T[tt]; iT2=iT*iT;

    for i in range(eqSteps):         # equilibrate
        mcmove(config, iT, D_data)           # Monte Carlo moves

    for i in range(mcSteps):
        mcmove(config, iT, D_data)
        Ene = calcEnergy(config)     # calculate the energy
        Mag = calcMag(config)        # calculate the magnetisation

        E1 = E1 + Ene
        M1 = M1 + Mag
        M2 = M2 + Mag*Mag
        E2 = E2 + Ene*Ene

    E_1[tt] = n1*E1
    M_1[tt] = n1*M1
    C_1[tt] = (n1*E2 - n2*E1*E1)*iT2
    X_1[tt] = (n1*M2 - n2*M1*M1)*iT

comm.send(E_1,dest=0,tag=rank)
comm.send(M_1,dest=0,tag=rank)
comm.send(C_1,dest=0,tag=rank)
comm.send(X_1,dest=0,tag=rank)
#print(rank,E_1)

if rank == 0:
    for i in range(0,size):
        E += comm.recv(source=i,tag=i)
        M += comm.recv(source=i,tag=i)
        C += comm.recv(source=i,tag=i)
        X += comm.recv(source=i,tag=i)
    with open('Energy.txt', 'w') as f:
	    for i in range(nt):
		print("{0:4d} {1:5f}".format(i,E[i]),file=f)
    with open('Polarization.txt', 'w') as f:
            for i in range(nt):
                print("{0:4d} {1:5f}".format(i,M[i]),file=f)
    with open('Specific_Heat.txt', 'w') as f:
            for i in range(nt):
                print("{0:4d} {1:5f}".format(i,C[i]),file=f)
    with open('Susceptibility.txt', 'w') as f:
            for i in range(nt):
                print("{0:4d} {1:5f}".format(i,X[i]),file=f)

    f = plt.figure(figsize=(18, 10)); # plot the calculated values

    sp =  f.add_subplot(2, 2, 1 );
    plt.scatter(T, E, s=50, marker='o', color='IndianRed')
    plt.xlabel("Temperature (T)", fontsize=20);
    plt.ylabel("Energy ", fontsize=20);         plt.axis('tight');

    sp =  f.add_subplot(2, 2, 2 );
    plt.scatter(T, abs(M), s=50, marker='o', color='RoyalBlue')
    plt.xlabel("Temperature (T)", fontsize=20);
    plt.ylabel("Magnetization ", fontsize=20);   plt.axis('tight');

    sp =  f.add_subplot(2, 2, 3 );
    plt.scatter(T, C, s=50, marker='o', color='IndianRed')
    plt.xlabel("Temperature (T)", fontsize=20);
    plt.ylabel("Specific Heat ", fontsize=20);   plt.axis('tight');

    sp =  f.add_subplot(2, 2, 4 );
    plt.scatter(T, X, s=50, marker='o', color='RoyalBlue')
    plt.xlabel("Temperature (T)", fontsize=20);
    plt.ylabel("Susceptibility", fontsize=20);   plt.axis('tight');
    #plt.show()
    plt.savefig("MC.png")
