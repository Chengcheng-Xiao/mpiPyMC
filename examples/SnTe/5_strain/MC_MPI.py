from __future__ import print_function
from __future__ import division
import numpy as np
from numpy.random import rand
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

#----------------------------------------------------------------------
##  BLOCK OF FUNCTIONS USED IN THE MAIN CODE
#----------------------------------------------------------------------
def initialstate(N,ps):
    ''' generates a random spin configuration for initial condition'''
    state = ps*np.random.randint(1,2, size=(N,N)) #1.51*(2*np.random.randint(2, size=(N,N))-1) #np.random.uniform(-2.98,2.98, size=(N,N))
    return state


def mcmove(config, beta, ps, A_data, B_data, C_data, D_data):
    '''Monte Carlo move using Metropolis algorithm '''
    for i in range(N):
        for j in range(N):
                a = np.random.randint(0, N)
                b = np.random.randint(0, N)
                # Nearest neighbour mean value
                nb = (config[(a+1)%N,b] + config[a,(b+1)%N] + config[(a-1)%N,b] + config[a,(b-1)%N])/4.
                # site polarization before change
                s =  config[a, b]
                # site polarization after change
                s1 = np.random.uniform(-2*ps,2*ps)
                # unitcell energy before and after
                cost_unitcell_orig = A_data/2.*s**2.+B_data/4.*s**4.+C_data/6.*s**6.
                cost_unitcell_change = A_data/2.*s1**2.+B_data/4.*s1**4.+C_data/6.*s1**6.
                # coupling energy before and after
                cost_coupling_orig = 4.*(D_data/2.)*(s-nb)**2.
                cost_coupling_change = 4.*(D_data/2.)*(s1-nb)**2.
		        # Calculate total energy change
                cost = (cost_unitcell_change + cost_coupling_change) - (cost_unitcell_orig + cost_coupling_orig)
                if cost < 0:
                    s = s1
                elif rand() < np.exp(-cost*beta):
                    s = s1
                config[a, b] = s
    return config

def calcMag(config):
    '''Magnetization of a given configuration'''
    mag = np.abs(np.sum(config))
    return mag

## change these parameters for a smaller (faster) simulation
nt      = 18         #  number of temperature points
N       = 30         #  size of the lattice, N x N
eqSteps = 8000       #  number of MC sweeps for equilibration
mcSteps = 4000       #  number of MC sweeps for calculation
A_data, B_data, C_data, D_data = -30.108, -5.928, 1.254, 8.642
ps         =  2.765

T       = np.linspace(800,1200, nt);
E,M,C,X = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
E_1,M_1,C_1,X_1 = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
n1, n2  = 1.0/(mcSteps*N*N), 1.0/(mcSteps*mcSteps*N*N)
# divide by number of samples, and by system size to get intensive values

#----------------------------------------------------------------------
#  MAIN PART OF THE CODE
#----------------------------------------------------------------------

for tt in range(int(nt*rank/size),int(nt*(rank+1)/size)):
    E1 = M1 = E2 = M2 = 0
    config = initialstate(N, ps)
    iT=1.0/(T[tt]*0.08617)

    for i in range(eqSteps):         # equilibrate
        mcmove(config, iT, ps, A_data, B_data, C_data, D_data)           # Monte Carlo moves

    for i in range(mcSteps):
        mcmove(config, iT, ps, A_data, B_data, C_data, D_data)
        Mag = calcMag(config)        # calculate the magnetisation

        M1 = M1 + Mag
    #average Polarization
    M_1[tt] = n1*M1

comm.send(M_1,dest=0,tag=rank)

if rank == 0:
    for i in range(0,size):
        # collect all data from MPI process
        M += comm.recv(source=i,tag=i)
    # Write data to Polarization.txt
    with open('Polarization.txt', 'w') as f:
            for i in range(nt):
                print("{0:4d} {1:5f}".format(i,M[i]),file=f)

# Plotting ....
    f = plt.figure(figsize=(18, 10)); # plot the calculated values

    plt.scatter(T, abs(M), s=50, marker='o', color='RoyalBlue')
    plt.xlabel("Temperature (T)", fontsize=20);
    plt.ylabel("Magnetization ", fontsize=20);   plt.axis('tight');

    plt.savefig("MC.png")

