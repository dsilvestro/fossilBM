# Scripts used in:
# Evolutionary history of New World monkeys revealed by molecular and fossil data
# by D Silvestro, M F Tejedor, M L Serrano-Serrano, O Loiseau, V Rossier, J Rolland, A Zizka, A Antonelli, N Salamin
# http://www.biorxiv.org/content/early/2017/08/18/178111
# doi: https://doi.org/10.1101/178111


# RUN ANALYSIS ON SIMULATED DATA and PLOT RESULTS
# --n 100000 sets the number of MCMC iterations to 100K (a higher number is likely necessary for convergence)
RScript mcmcFossilBM.R --wd /path_to_mcmcFossilBM_file --n 100000
# Additional commands for simulations
"""
flag    |  type     |  default |  help
--t     |  integer  |  50      |  Tree size 
--v     |  double   |  16      |  Magnitude rate shift
--r     |  double   |  0       |  Baseline rate
--s     |  integer  |  0       |  n. shifts sig2 
--mu    |  double   |  NA      |  Baseline trend parameter mu0 
--sm    |  integer  |  0       |  n. shifts mu0
--rda   |  integer  |  0       |  if 1 use existing RDA file 
--qfos  |  integer  |  20      |  Mean number of simulated fossils
--j     |  integer  |  1       |  Simulation number (goes in output file names) 
"""

# RUN EMPIRICAL DATA
# latitude, values rescaled (x 0.05) to improve convergence
RScript mcmcFossilBM.R --wd /path_to_tree_and_data_files --n 100000 --tfile platyrrhineFBD.trees --dfile platyrrhine_latitude.txt --rescale 0.05 --out LAT 
# body mass, values log10 transformed
RScript mcmcFossilBM.R --wd /path_to_tree_and_data_files --n 100000 --tfile platyrrhineFBD.trees --dfile platyrrhine_bodymass.txt --log 10 --out BM


