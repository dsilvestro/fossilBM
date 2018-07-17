## FossilBM - Bayesian inference of trait evolution under a “fossilized” Brownian motion (FBM)

FossilBM uses phylogenies of extant and extinct species (typically obtained through tip-dating) and quantitative traits to infer the mode of evolution and ancestral states at all nodes in the tree.
The FBM model includes rate parameters describing how fast a trait evolves and trend parameters determining whether evolution is directional. When the trend parameter equals 0, the model reduces to a standard Brownian model of neutral evolution. Both the rate and trend parameters can very across the branches in the tree.

The algorithm simultaneously estimates:  
1. the rates of phenotypic evolution and their variation across clades  
2. the trends of phenotypic evolution and their variation across clades  
3. the ancestral states at all internal nodes


The methods implemented here are described in:
Early arrival and climatically-linked geographic expansion of New World monkeys from tiny African ancestors
by D Silvestro, M F Tejedor, M L Serrano-Serrano, O Loiseau, V Rossier, J Rolland, A Zizka, A Antonelli, N Salamin ([Systematic Biology, 2018](https://doi.org/10.1093/sysbio/syy046)).


### Run analysis on simulated data and plot results

 `RScript mcmcFossilBM.R --wd /path_to_mcmcFossilBM_file --n 100000`
 
where `--n 100000` sets the number of MCMC iterations to 100K (a higher number is likely necessary for convergence)

#### Additional commands for simulations:  

```
flag    |  type     |  default |  help
--t     |  integer  |  50      |  Tree size 
--v     |  double   |  16      |  Magnitude rate shift
--r     |  double   |  0       |  Baseline rate parameter (sig2)
--s     |  integer  |  0       |  n. shifts sig2 
--mu    |  double   |  NA      |  Baseline trend parameter (mu0) 
--sm    |  integer  |  0       |  n. shifts mu0
--rda   |  integer  |  0       |  if 1 use existing RDA file 
--qfos  |  integer  |  20      |  Mean number of simulated fossils
--j     |  integer  |  1       |  Simulation number (goes in output file names) 
```
___
### Run the analysis on empirical data
This command analyzes mid latitude of species ranges as a trait after specifying the tree file (`--tfile`; see example files in the repository) and a trait table (`--dfile`):

`RScript mcmcFossilBM.R --wd /path_to_tree_and_data_files --n 100000 --tfile platyrrhineFBD.trees --dfile platyrrhine_latitude.txt --rescale 0.05 --out LAT `

The latitude values are rescaled (x0.05) using the command `--rescale` to improve the convergence of the MCMC. The command `--out LAT` adds "LAT" to the output file names.

To analyze body mass data we need to log transform the data using the command `--log 10`:
`RScript mcmcFossilBM.R --wd /path_to_tree_and_data_files --n 100000 --tfile platyrrhineFBD.trees --dfile platyrrhine_bodymass.txt --log 10 --out BM`


