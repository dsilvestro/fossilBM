## FossilBM - Bayesian inference of trait evolution under a “fossilized” Brownian motion (FBM)

FossilBM uses phylogenies of extant and extinct species (typically obtained through tip-dating) and quantitative traits to infer the mode of evolution and ancestral states at all nodes in the tree.
The FBM model includes rate parameters describing how fast a trait evolves and trend parameters determining whether evolution is directional. When the trend parameter equals 0, the model reduces to a standard Brownian model of neutral evolution. Both the rate and trend parameters can very across the branches in the tree.

The algorithm simultaneously estimates:  
1. the rates of phenotypic evolution and their variation across clades  
2. the trends of phenotypic evolution and their variation across clades  
3. the ancestral states at all internal nodes

An example of a fossilBM setup and analysis is available in the file [`run_FBM.R`](https://github.com/dsilvestro/fossilBM/blob/master/run_FBM.R). 


The methods implemented here are described in [Silvestro et al. 2018 Systematic Biology](https://doi.org/10.1093/sysbio/syy046).
The BM model with time-varying trend is described in [Zhang et al. 2022 Systematic Biology](https://doi.org/10.1093/sysbio/syab065).  
The BM model with constrained trend paramenbters is described in [Farina et al. 2023 Proc Roy Soc B](https://doi.org/10.1098/rspb.2023.1099). 
  

#### NOTE
The following libraries are required to run FossilBM (`fossilBM_lib.R`): scales, phytools, geiger, adephylo, plotrix. 
Older versions of fossilBM model are available [here](https://github.com/dsilvestro/fossilBM/releases/tag/v0.2). 

