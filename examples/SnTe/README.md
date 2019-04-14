## Discrepancies

1. Using mpiPyMC, I am able to plot the phase transition for SnTe as published in [paper](https://aip.scitation.org/doi/10.1063/1.4996171). The T_c I got is `160K @ 0%_strain` and `930K @ 5%_strain`. The unstrained one is identical to the published result, the 5% strained one is deviated a little to the published result (1023K). The cause of this deviation may be the convergence criterial, which is typically higher in high temperature region.
