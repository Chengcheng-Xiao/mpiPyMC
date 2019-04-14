# mpiPyMC

A python based, MPI enabled, Monte-Carlo calculation of  2D ferroelectric system using Metropolis algorithm.

Before getting into things, you might want to check out these papers:

1.  [Phys. Rev. Lett. 117, 097601 (2016)](https://link.aps.org/doi/10.1103/PhysRevLett.117.097601)
2.  [Phys. Rev. B 97, 144104 (2018)](https://link.aps.org/doi/10.1103/PhysRevB.97.144104)
3.  [Appl. Phys. Lett. 111, 132904 (2017)](https://aip.scitation.org/doi/10.1063/1.4996171)

## Getting Started

These instructions will get you a copy of the project up and running on your machine.

### Prerequisites

For the script to work, you need to have an valid installation of `python` (2.7.x or 3.x both work).
Also, `numpy`, `matplotlib` and `mpi4py` package are needed, you can install them by pip:
```
pip install matplotlib numpy mpi4py
```
or by conda
```
conda install matplotlib numpy mpi4py
```
if you use a supercomputer and don't have enough privilege:

1. install anaconda by downloading from [here](https://www.anaconda.com/download/) and upload it to your directory.
2. using queue system to install anaconda by `chmod 755 anaconda*.sh && ./anaconda*.sh`
3. install `numpy`, `h5py` and `mpi4py` by download them from [here](https://anaconda.org/anaconda/repo) upload them as well.
4. manually install package by `conda install *package_name*.tar.bz2`

### Installing

Python script does not need manual compilation, so the installation is very easy.

1. Download the script:
```
wget https://github.com/Chengcheng-Xiao/mpiPyMC/blob/master/MC_MPI.py
```

2. give correct permission:
```
chmod 755 MC_MPI.py
```

3. change parameters inside the code and run it by:
```
./MC_MPI.py
```

## Adjustable variables

All adjustable variables are currently located inside the code. This will change in the future revisions.

| NAME                   | REQURE                                     |
|:----------------------:|:------------------------------------------:|
| `nt`                   | [number of Temperature point]                  |
| `N`                    | [number of cell in onedirection, total number of cell  =N*N]              |
| `eqSteps`              | [number of MC steps to get to equilibrium]              |
| `mcSteps`              | [number of MC steps to use for average]                        |
| `A_data``B_data``C_data``D_data`| [model data]                    |
| `ps`                   | [spontaneous polarization value]             |
| `T`                    | [temperature range]             |


## Example
In the `example` folder, I have included three examples for mentioned reference papers. The systems are: `SnSe`, `b-GeSe` and `SnTe`.

I was not able to reproduce the result for `SnSe` and `b-GeSe`.
The calculated `SnTe` results agree with [said paper](https://aip.scitation.org/doi/10.1063/1.4996171).

Analysis of possible errors are presented inside each folder separately.

## Note
This code is based on [🔗LINK](https://rajeshrinet.github.io/blog/2014/ising-model/)

## Future development plane
1. Add total energy plot, specific heat plot and susceptibility plot function.
2. Better MPI implementation?
3. Snap shot of specific one (or several) MC steps.

## License
  This project is licensed under the GNU License - see the `LICENSE.md` for details
