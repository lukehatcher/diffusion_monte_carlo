# Diffusion Monte Carlo
Diffusion Monte Carlo (DMC) is an algorithm commonly used in theoretical chemistry to calculate the zero point energies of a molecular system. As the name implies, the algorithm utilizes Monte Carlo random sampling. Common systems studied with this algorithm are small sized water clusters ((H2O)N). A method called *descendant weighting* can be concurrently woven into the DMC algorithm in order to simultaneously obtain the wave functions of a system. This repository includes a 1-dimensional and an N-dimensional version of the algorithm, both implementing descendant weighting, under dmc_1D_DW.py and dmc_ND_DW.py respectively. These implementations are done using a functional approach and would benefit from a more generalized, OOP approach.

## Mathematical motivation
Please see `dmc_derivation_and_background.pdf` for a PDF or `dmc_derivation_and_background.tex` for a LaTeX version. Typo corrections to come. 
