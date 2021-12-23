# mpiPyMC

A python based, MPI enabled, Monte-Carlo calculation of  2D Ising system using Metropolis algorithm.

Before getting into things, you might want to check out these papers:

1.  [Phys. Rev. Lett. 117, 097601 (2016)](https://link.aps.org/doi/10.1103/PhysRevLett.117.097601)
2.  [Phys. Rev. B 97, 144104 (2018)](https://link.aps.org/doi/10.1103/PhysRevB.97.144104)
3.  [Appl. Phys. Lett. 111, 132904 (2017)](https://aip.scitation.org/doi/10.1063/1.4996171)

## Getting Started

These instructions will get you a copy of the project up and running on your machine.

### Prerequisites

For the script to work, you need to have an valid installation of `python` (2.7.x or 3.x both work), and a MPI installation:
```
openmpi:https://www.open-mpi.org/
MPICH:https://www.mpich.org/
intelmpi:https://software.intel.com/en-us/mpi-library
```

Also, `numpy`, `matplotlib`, `tqdm` and `mpi4py` package are needed, you can install them by pip:
```
pip install matplotlib numpy mpi4py tqdm
```
or by conda
```
conda install matplotlib numpy mpi4py tqdm
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
wget https://github.com/Chengcheng-Xiao/mpiPyMC/blob/master/MC_MPI*.py
```

2. give correct permission:
```
chmod 755 MC_MPI*.py
```

3. change parameters inside the code and run it by:
```
mpirun -np XX MC_MPI*.py
```
`XX` is for number of processes (please keep XX < the number of temperature points).

### Adjustable spin (polarization vector)
This project provide two different script for two different tasks:
1. `MC_MPI.py`: for adjustable spin. Utilize Sigma^4 Model:
```
E=∑[(A/2)*P_i^2+(B/4)*P_i^4+(C/6)*P_i^6]+∑[(D/2)*(P_i-<P_j>)^2]
```
where i and j are nearest neighbor.

2. `MC_MPI_Ising.py`: for spin = +1 or -1. Utilize Heisenberg Model:
```
E=∑[(D)*P_i*P_j]
```
where i and j are nearest neighbor.


## Input File

All adjustable variables should be written in `input.MC` file.

FOR `MC_MPI.py`:

| NAME                   |  DEFAULT     |MEANING                                     |
|:----------------------:|:---------------:|:---------------------------:|
| `nt`                   |        18       |[number of Temperature point, should be integer number of the CPU used]                  |
| `N`                    |        16       |[number of cell in onedirection, total number of cell  =N*N]              |
| `eqSteps`              |       2000      |[number of MC steps to get to equilibrium]              |
| `mcSteps`              |       2000      |[number of MC steps to use for average]                        |
| `A_data`               |      -8.021     |[model data]                    |
| `B_data`               |      0.620      |[model data]                    |
| `C_data`               |      0.441      |[model data]                    |
| `D_data`               |        4.5      |[model data]                    |
| `ps`                   |     1.867       |[spontaneous polarization value]             |
| `T_low`                |        0.01     |[temperature range]             |
| `T_high`               |        300      |[temperature range]             |


FOR `MC_MPI_Ising.py`:

| NAME                   |  DEFAULT     |MEANING                                     |
|:----------------------:|:---------------:|:---------------------------:|
| `nt`                   |        18       |[number of Temperature point, should be integer number of the CPU used]                  |
| `N`                    |        16       |[number of cell in onedirection, total number of cell  =N*N]              |
| `eqSteps`              |       8000      |[number of MC steps to get to equilibrium]              |
| `mcSteps`              |       4000      |[number of MC steps to use for average]                        |
| `D_data`               |        1.0      |[model data]                    |
| `T_low`                |        1.53     |[temperature range]  assuming K_b=1           |
| `T_high`               |        3.28      |[temperature range] assuming K_b=1            |



### Output
`MC_MPI.py` will output picture `MC.png` and data `Polarization.txt`, containing Polarization vs Temperature plot and the raw data to generate the plot, respectively.

`MC_MPI_Ising.py` will output picture `MC.png` containing total energy plot, specific heat plot, susceptibility plot and Polarization vs Temperature plot.
as well as their raw data `Energy.txt`, `Polarization.txt`, `Specific_Heat.txt`, `Susceptibility.txt`

## Example
### `MC_MPI` example
In the `example` folder, I have included three examples for mentioned reference papers. The systems are: `SnSe`, `ß-GeSe`, `SnTe` and `GeTe`.

I was not able to reproduce the result for `SnSe` and `ß-GeSe`.
The calculated results for`SnTe` and `GeTe` agree with [said paper](https://aip.scitation.org/doi/10.1063/1.4996171).

Analysis of possible errors are presented inside each folder separately. In short, the definition of nearest neighbor is very important. However, no explanations were made in these referenced paper. The coupling coefficient D should be calculated via changing one cell's polarization in a supercell configuration. The coupling plot should be calculated by subtracting the site energy from the total energy.

### `MC_MPI_Ising` example
In the `example_Ising` folder, I have included three calculations of Ising model. each with different cell size: 16X16, 30X30 and 50X50.

## Note
This code is based on rajeshrinet's work: [🔗LINK](https://rajeshrinet.github.io/blog/2014/ising-model/)

## Future development plan
1. Add total energy plot, specific heat plot and susceptibility plot function for `MC_MPI.py`.
2. Better MPI implementation with respect to CPU number, or use OpenMP?
3. Snap shot of specified one (or several) MC step(s).
4. Hysteresis.

## License
  This project is licensed under the MIT License - see the [LICENSE.md](./LICENSE.md) for details
