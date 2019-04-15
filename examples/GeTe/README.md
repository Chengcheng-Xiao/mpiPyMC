## Discrepancies

1. Using mpiPyMC, I am able to plot the phase transition for GeTe as published in [paper](https://aip.scitation.org/doi/10.1063/1.4996171). The T_c I got is `570K @ 0%_strain` and `1100K @ 5%_strain`. The unstrained one is identical to the published result, the 5% strained one is deviated a little to the published result (1350K). The cause of this deviation may be the convergence criterial, which is typically higher in high temperature region and require larger supercell.

## Parameters for `0_strain`
| NAME                   | DATA                                     |
|:----------------------:|:------------------------------------------:|
| `nt`                   | 18                 |
| `N`                    | 30              |
| `eqSteps`              | 8000              |
| `mcSteps`              | 4000                        |
| `A_data``B_data``C_data``D_data`| -14.482, -2.155, 0.315, 3.364                    |
| `ps`                   | 3.276             |
| `T`                    | 400, 800             |

## Parameters for `5_strain`
| NAME                   | DATA                                     |
|:----------------------:|:------------------------------------------:|
| `nt`                   | 18                 |
| `N`                    | 30              |
| `eqSteps`              | 8000              |
| `mcSteps`              | 4000                        |
| `A_data``B_data``C_data``D_data`| -24.997, -9.493, 0.920, 4.845                    |
| `ps`                   | 3.493             |
| `T`                    | 900,1300             |
