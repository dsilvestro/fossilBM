## FossilBM - Bayesian inference of trait evolution under a “fossilized” Brownian motion (FBM)

FossilBM uses phylogenies of extant and extinct species (typically obtained through tip-dating) and quantitative traits to infer the mode of evolution and ancestral states at all nodes in the tree.
The FBM model includes rate parameters describing how fast a trait evolves and trend parameters determining whether evolution is directional. When the trend parameter equals 0, the model reduces to a standard Brownian model of neutral evolution. Both the rate and trend parameters can very across the branches in the tree.

The algorithm simultaneously estimates:  
1. the rates of phenotypic evolution and their variation across clades  
2. the trends of phenotypic evolution and their variation across clades  
3. the ancestral states at all internal nodes


The methods implemented here are described in:
*Early arrival and climatically-linked geographic expansion of New World monkeys from tiny African ancestors* by D. Silvestro, M. F. Tejedor, et al. ([Systematic Biology, 2018](https://doi.org/10.1093/sysbio/syy046)).
  
   
___
### Run the analysis on empirical data
To run a fossiBM analysis you should provide at least three arguments:   

1) a tree file in nexus format  
2) a simple text file including taxa names and trait values as two tab-separated columns  
3) working directory (where the input files are and where the output files will be saved).  

To run an analysis run in a Terminal:
 `RScript mcmcFossilBM.R --wd /path_to_files --tfile tree_file --dfile trait_file`


#### Command list:  

```
flag        |  type     |  default |  help
--tfile     |  string   |  none      |  Tree file 
--dfile     |  string   |  none      |  Trait file
--wd        |  string   |  none      |  Working directory
--n         |  integer  |  250000    |  n. MCMC iterations
--tindex    |  integer  |  1         |  index of tree analyzed (when tree file includes multiple trees)
--out       |  string   |  none      |  add suffix to output files
--rescale   |  integer  |  1         |  multiply trait values to rescale them (can help convergence)
--log       |  integer  |  0         | log transform trait values - options: 10 (log10), 1 (ln), 0 (no transformation)
--rmF       |  integer  |  0         |  If set to 1 removes all fossil tips from the tree
--useBMT    |  integer  |  1         |  If set to 0 runs BM model without trend
--constRate |  integer  |  0         |  If set to 1 runs BMT model with constant rate
--pfile     |  string   | none       |  Specify partition file for clade-specific rates and trends
```  


The following commands replicate the analyses presented in [Silvestro et al. 2018 Syst Biol](https://doi.org/10.1093/sysbio/syy046) (input files are provided in the [GitHub repository](https://github.com/dsilvestro/fossilBM)).

This command analyzes mid latitude of species ranges as a trait:

`RScript mcmcFossilBM.R --wd /path_to_tree_and_data_files --n 100000 --tfile platyrrhineFBD.trees --dfile platyrrhine_latitude.txt --rescale 0.05 --out LAT `

The latitude values are rescaled (x0.05) using the command `--rescale` to improve the convergence of the MCMC. The command `--out LAT` adds "LAT" to the output file names.

To analyze body mass data we need to log transform the data using the command `--log 10`:
`RScript mcmcFossilBM.R --wd /path_to_tree_and_data_files --n 100000 --tfile platyrrhineFBD.trees --dfile platyrrhine_bodymass.txt --log 10 --out BM`


The results are **plotted automatically** at the end of the run, and the plots are saved in a PDF file. 

#### Running under models with fixed tree partitions ####

You can constrain the analysis to use pre-define partitions, i.e. to estimate independent rates and trends in different clades of the tree. The clades are defined in a text file (see 'clade_partitions.txt' example) which specify two taxa identifying the clade(s) of interest through their most recent common ancestor, e.g.

```
Alouatta_belzebul	      Stirtonia_victoriae
Mazzonicebus_almendrae    Cacajao_calvus
```

The command `--pfile clade_partitions.txt` is then used to enforce these partition. Hence, the command will look like: 

`RScript mcmcFossilBM.R --tfile platyrrhine_FBD.trees --dfile platyrrhine_bodymass.txt --log 10 --wd /path_to_tree_and_data_files --pfile clade_partitions.txt `

Note that enforcing partitions will automatically switch off the BDMCMC option to estimate the number of rate and trend shifts. 






---
### Additional commands for simulations
The following commands can be used to run analyses on simulated data.

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