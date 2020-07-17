# Diffusion Monte Carlo
Diffusion Monte Carlo (DMC) is an algorithm commonly used in theoretical chemistry to calculate the zero point energies of a molecular system. As the name implys, the algorithm hinges on Monte Carloesque random sampling. Common systems studied with this algorithm are small sized water clusters ((H<sup>2</sup>O)<sup>N</sup>). A method called descendant weighting can be concurrently woven into the DMC algorithm in order to simultaneously obtain the wave functions of system of study. This repository includes a 1-dimensional and an N-dimensional version of the algorithm, both with implemented descendant weighting, under `dmc_1D_DW.py` and `dmc_ND_DW.py` respectively. These implementations are done using funcitonal programming, and would benefit from a more generalized, oop approach. 

## Mathematical motivation
Please compile `dmc_derivation_and_background.tex` for a lengthy overview. PDF with much needed typo corrections coming soon. 
