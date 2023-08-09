
Code related to the paper "Bias reduction and robustness for gaussian longitudinal data analysis" (E. Bortolato, E. C. Kenne Pagui, 2023)
to appear in Journal Of Statistical Computation And Simulation.

# bren (Bias Reduction in the Equicorrelated multivariate Normal)
The function "bren" fits regression for gaussian equicorrelated data using adjusted score functions (bias correction, mean bias reduction, median bias reduction, Jeffreys (power) penalization.
Jeffreys (power) penalization method is robust (i.e. produces reliable confidence intervals) to misspecified correlation structures.

# Stroke example
Example for using "bren" with the Stroke dataset (Dobson e Barnett, 2008), available in the \texttt{R} package \text{MLGdata} (Sartori et al.) on CRAN.
This file also contains the code for reproducing the simulation study using this dataset in the paper

# simulation study Scenarios 1 and 2
Code for reproducing the simulation study in section 4 and producing tables and figures
