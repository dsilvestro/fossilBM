## FossilBM - Bayesian inference of trait evolution under a “fossilized” Brownian motion (FBM)

FossilBM uses phylogenies of extant and extinct species (typically obtained through tip-dating) and quantitative traits to infer the mode of evolution and ancestral states at all nodes in the tree.
The FBM model includes rate parameters describing how fast a trait evolves and trend parameters determining whether evolution is directional. When the trend parameter equals 0, the model reduces to a standard Brownian model of neutral evolution. Both the rate and trend parameters can very across the branches in the tree.

The algorithm simultaneously estimates:  
1. the rates of phenotypic evolution and their variation across clades  
2. the trends of phenotypic evolution and their variation across clades  
3. the ancestral states at all internal nodes


The methods implemented here are described in Silvestro, Tejedor, et al. ([Systematic Biology, 2018](https://doi.org/10.1093/sysbio/syy046)).
  

#### NOTE
Two versions of the fossilBM model are available: version.1 is launched from a terminal window using `RScript mcmcFossilBM.R [additional commands]`. A newer implementation loads all the necessary functions from the file `fossilBM_lib.R` and can be run more interactively from an R console. Current and future developments of fossilBM are implemented in `fossilBM_lib.R` only. 


The following libraries are required to run FossilBM (`fossilBM_lib.R`): scales, phytools, geiger, adephylo. An example of a fossilBM setup is available in the file `run_FBM.R`. 
