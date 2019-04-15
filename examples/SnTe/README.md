## Discrepancies

1. Using mpiPyMC, I am able to plot the phase transition for SnTe as published in [paper](https://aip.scitation.org/doi/10.1063/1.4996171). The T_c I got is `160K @ 0%_strain` and `930K @ 5%_strain`. The unstrained one is identical to the published result, the 5% strained one is deviated a little to the published result (1023K). The cause of this deviation may be the convergence criterial, which is typically higher in high temperature region and require larger supercell.

## Parameters for `0_strain`
| NAME                   | DATA                                     |
|:----------------------:|:------------------------------------------:|
| `nt`                   | 18                 |
| `N`                    | 30              |
| `eqSteps`              | 8000              |
| `mcSteps`              | 4000                        |
| `A_data``B_data``C_data``D_data`| -8.021, 0.620, 0.441, 4.5                    |
| `ps`                   | 1.867             |
| `T`                    | 0.01, 300             |

## Parameters for `5_strain`
| NAME                   | DATA                                     |
|:----------------------:|:------------------------------------------:|
| `nt`                   | 18                 |
| `N`                    | 30              |
| `eqSteps`              | 80000              |
| `mcSteps`              | 100000                        |
| `A_data``B_data``C_data``D_data`| -30.108, -5.928, 1.254, 8.642                    |
| `ps`                   | 2.765             |
| `T`                    | 800,1200             |
