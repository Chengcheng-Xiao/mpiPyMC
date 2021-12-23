#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from tqdm.auto import tqdm
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
def read_input():

    ''' a subroutine to get AO and MO number and number of kpoints and otehr information from input.woops'''

    dataset={}
    file = open('input.MC', "r")
    #Default
    data = file.readlines()
    nt      = 18         #  number of temperature points
    N       = 16         #  size of the lattice, N x N
    eqSteps = 2000       #  number of MC sweeps for equilibration
    mcSteps = 2000       #  number of MC sweeps for calculation
    A_data, B_data, C_data, D_data = -15.261, 1.8389, 1.7137, 8.0904
    ps         =  1.59
    T_low      =  0.01
    T_high     =  400

    for line in data:
        key, value = line.split("=")
        dataset[key.strip()] = value.strip()
    # Read data
    nt          = int(dataset["nt"])
    N           = int(dataset["N"])
    eqSteps     = int(dataset["eqSteps"])
    mcSteps     = int(dataset["mcSteps"])
    A_data      = float(dataset["A_data"])
    B_data      = float(dataset["B_data"])
    C_data      = float(dataset["C_data"])
    D_data      = float(dataset["D_data"])
    ps          = float(dataset["ps"])
    T_low       = float(dataset["T_low"])
    T_high      = float(dataset["T_high"])
    return nt, N, eqSteps, mcSteps, A_data, B_data, C_data, D_data, ps, T_low, T_high

def initialstate(N,ps):
    ''' generates a random spin configuration for initial condition'''
    state = ps*np.random.randint(1,2, size=(N,N)) #1.51*(2*np.random.randint(2, size=(N,N))-1) #np.random.uniform(-2.98,2.98, size=(N,N))
    return state


def mcmove(config, beta, ps, A_data, B_data, C_data, D_data):
    '''Monte Carlo move using Metropolis algorithm '''
    for i in range(N*N):
        a = np.random.randint(0, N)
        b = np.random.randint(0, N)
        # Nearest neighbour mean value
        nb = (config[(a+1)%N,b] + config[a,(b+1)%N] + config[(a-1)%N,b] + config[a,(b-1)%N])/4.
        # site polarization before change
        s =  config[a, b]
        # site polarization after change, note here I constrain the max polarization to be 2*ps.
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

#########################################################################
# Here comes the Model parameters
#########################################################################
nt, N, eqSteps, mcSteps, A_data, B_data, C_data, D_data, ps, T_low, T_high = read_input()
#########################################################################
T               = np.linspace(T_low, T_high, nt);
E,M,C,X         = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
E_1,M_1,C_1,X_1 = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
n1, n2          = 1.0/(mcSteps*N*N), 1.0/(mcSteps*mcSteps*N*N)
# divide by number of samples, and by system size to get intensive values

#----------------------------------------------------------------------
#  MAIN PART OF THE CODE
#----------------------------------------------------------------------

# exit if we have more processes then the number of temperature points
if nt < size-1:
    if rank == 0:
        print('**ERROR: More processes then the number of temperatrue points.')
    exit()

# automatically convert number of temperature steps to be devideable by "size"
if rank == 0:
    print("                  _ ____        __  __  ____ \n"
    "  _ __ ___  _ __ (_)  _ \ _   _|  \/  |/ ___|\n"
    " | '_ ` _ \| '_ \| | |_) | | | | |\/| | |    \n"
    " | | | | | | |_) | |  __/| |_| | |  | | |___ \n"
    " |_| |_| |_| .__/|_|_|    \__, |_|  |_|\____|\n"
    "           |_|            |___/              \n")
    pbar = tqdm(total=nt*(eqSteps+mcSteps),desc="total_% ")
    # for tt in tqdm(range(int(nt*rank/size),int(nt*(rank+1)/size)),desc="total_% ",ncols=100):
    #     E1 = M1 = E2 = M2 = 0
    #     config = initialstate(N, ps)
    #     iT=1.0/(T[tt]*0.08617)
    #
    #     for i in range(eqSteps):
    #     # for i in tqdm(range(eqSteps),desc="Eq_steps",leave=False,ncols=100):         # equilibrate
    #         mcmove(config, iT, ps, A_data, B_data, C_data, D_data)           # Monte Carlo moves
    #
    #     for i in range(mcSteps):
    #     # for i in tqdm(range(mcSteps),desc="MC_steps",leave=False,ncols=100):
    #         mcmove(config, iT, ps, A_data, B_data, C_data, D_data)
    #         Mag = calcMag(config)        # calculate the magnetisation
    #
    #         M1 = M1 + Mag
    #     #average Polarization
    #     M_1[tt] = n1*M1
    #
    # comm.send(M_1,dest=0,tag=rank)

if rank != 0:
    for tt in range(int(nt*(rank-1)/(size-1)),int(nt*((rank-1)+1)/(size-1))):
        E1 = M1 = E2 = M2 = 0
        config = initialstate(N, ps)
        iT=1.0/(T[tt]*0.08617)

        for i in range(eqSteps):         # equilibrate
            mcmove(config, iT, ps, A_data, B_data, C_data, D_data)           # Monte Carlo moves
            comm.send(None,dest=0,tag=0)

        for i in range(mcSteps):
            mcmove(config, iT, ps, A_data, B_data, C_data, D_data)
            Mag = calcMag(config)        # calculate the magnetisation

            M1 = M1 + Mag
            comm.send(None,dest=0,tag=0)
        #average Polarization
        M_1[tt] = n1*M1

        # print(str(rank)+'fuck',flush=True)

    comm.send(M_1,dest=0,tag=rank)

if rank == 0:
    # for i in tqdm(range(1,size),desc="total_% sub ",ncols=100):
    remaining = size - 1
    while remaining > 0:
        s = MPI.Status()
        comm.Probe(status=s)
        if s.tag == 0:
            comm.recv(tag=0)
            pbar.update()
        else:
            # collect all data from MPI process
            M += comm.recv(source=s.tag,tag=s.tag)
            remaining -= 1
    pbar.close()
    # Write data to Polarization.txt
    with open('Polarization.txt', 'w') as f:
            for i in range(nt):
                print("{0:4d} {1:5f}".format(i,M[i]),file=f)

# Plotting ....
    f = plt.figure(figsize=(18, 10)); # plot the calculated values

    plt.scatter(T, abs(M), s=50, marker='o', color='RoyalBlue')
    plt.xlabel("Temperature (K)", fontsize=20);
    plt.ylabel("Magnetization ", fontsize=20);   plt.axis('tight');

    plt.savefig("MC.png")
    print("\n\n Done.")
