# FossilBM tutorial using trait states
_Bruna M. Farina_  
_2024-11-19_

 FossilBM can sample evolutionary parameters and ancestral states across phylogenies. It implements Brownian motion evolutionary models and allows evolutionary rates to vary across the phylogeny.

 Here is an example of how to setup a FossilBM analysis to infer body mass evolution of platyrrhine using state-dependent μ (trend parameter) and σ (rate parameter). Data used in this tutorial are available [here](https://github.com/dsilvestro/fossilBM/tree/master/example_files).

 Let's get started! :computer: 

## Loading and preparing your data
Let's start by loading and tidying up the necessary data.

````R
library(corHMM)

treefile  <- read.nexus("example_files/platyrrhine_FBD.trees")
tree <- treefile[[1]] # select only one tree 

datafile <- read.table("example_files/platyrrhine_bodymass.txt")
colnames(datafile) <- c("taxon_names", "body_mass") #giving the table meaningful column names

#tree and data have different number of species but they should match
datafile <- datafile[datafile$taxon_names %in% tree$tip.label,]
tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% datafile$taxon_names)])

rownames(datafile) <- datafile$taxon_names # assign taxon_names as rownames
datafile$taxon_names <- NULL
````

## Model testing and ancestral state estimates

For this example we are going to create a fake trait data.

````R
states <- round(runif(n=112, min=1, max=2), 0) #create a fake 2-state trait as an example
trait_data <- as.data.frame( cbind(datafile$taxon_name, states)) # create fake trait dataframe with taxon names
````

In order to model the evolution of a certain trait, we should fit different Markov Models. Here we are going to use [corHMM](https://cran.r-project.org/web/packages/corHMM/index.html) for that. Below we fit equal rates (ER) and all rates different (ARD) models (but the package has other models, including customized ones!)

````R
ER_model <- corHMM(tree, trait_data, model="ER", rate.cat=1, root.p=NULL)

ARD_model <- corHMM(tree, trait_data, model="ARD", rate.cat=1, root.p=NULL)
````
Once we fit the models, we can compare them using their AICc scores. Then we perform an ancestral state estimates using the transition rates estimated for the best model (lower AICc). Let's use the ARD model here:

````R
params <- sapply(na.omit(c(ARD_model$index.mat)),
                 function(x) na.omit(c(ARD_model$solution))[x])

lik_states <- ancRECON(phy= ARD_model$phy, data=ARD_model$data,
p = params, method = "joint", rate.mat = ARD_model$index.mat, rate.cat = ARD_model$rate.cat, root.p = NULL)
````

And let's create the trait table for running FossilBM.
````R
#get tip states
tip_states <- apply(ARD_model$tip.states, 1, which.max)
tip_state_tbl <- cbind(1:nrow(tip_trait), tip_states)

#anc states
anc_state_trait <- as.data.frame(cbind(seq(1:length(lik_states$lik.anc.states)), lik_states$lik.anc.states))
colnames(tip_state_tbl) = colnames(anc_state_trait)

# combine vectors of tip and anc state values
state_tbl <- rbind(tip_state_tbl, anc_state_trait[2:dim(anc_state_trait)[1],])
````

## FossilBM analysis
Now you should be able run FossilBM using the states prepared above!  
Let's create a FossilBM abject with the input data:

````R
fbm_obj <- read_and_transform_data(tree_obj=tree, 
data_obj=datafile, state_tbl=state_tbl, z_tranform=TRUE)
````
Arguments for ````read_and_transform_data```` function include:

|  Argument            |  What it does                                                                                                |
|---------------------:|-----------------------------------------------------------------------------------------------------------|
|_log_trait_data_      |if set to 0 trait is not log-transformed                                                                   |
|_drop_na_             |if set to TRUE NAs will be dropped otherwise they are imputed                                              |
|_rescale_trait_data_  |multiplier to rescale the trait data                                                                       |
|_root_calibration_    |normal prior (mean,sd) on root state                                                                       |
|_partition_file_      |independent rate and trends parameters for each partition (example [here](https://github.com/dsilvestro/fossilBM/tree/master/example_files))

And we can run the MCMC:

````R
output_file <- "fbm_mcmc.log"

res <- run_mcmc(fbm_obj, logfile=output_file,
                useTrend=TRUE, # BM model with trend
                linTrend=FALSE, per_branch_parameters = FALSE,
                log_anc_states = FALSE,
                ngen=10000, sample=250, print_freq=1000, 
                ) # mcmc settings
````

Other optional arguments for ````run_mcmc```` function include:

|  Argument     |  What it does                                                                                                          |
|--------------:|--------------------------------------------------------------------------------------------------------------------|
|_useTrend_     |set to FALSE to run a BM model with no trend                                                                        |
|_constRate_    |set to TRUE to run a BM model with constant rates (but variable trends)                                             |
|_linTrend_     |set to TRUE to run a BM model in which the trend paramater itself follows a linear trend                            |
|_verbose_      |set to TRUE to see rate and trend parameters                                                                        |
|_update_mu0_   |allows for controlling trend paramenter (e.g. update_mu0 = c(2, 3), three state trait, but keeping state 1 constant)|
