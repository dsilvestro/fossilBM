# define input files
setwd("/Users/dsilvestro/Software/fossilBM/example_files") # set working directory
treefile <- "pruned_tree.tre"
tindex <- 1 # specify which tree should be analyzed if the tree file includes multiple trees
datafile <- "platyrrhine_bodymass.txt"

# load FBM library
source("/Users/dsilvestro/Software/fossilBM/fossilBM_lib.R")

# process and match input data 
fbm_obj <- read_and_transform_data(treefile, datafile, log_trait_data=10, partition_file="nested_partitions.txt")



# run MCMC
output_file <- "fbm_mcmc.log"
run_mcmc(fbm_obj, logfile=output_file, ngen = 5000, sample = 50, print_freq=100, 
		useTrend = T, linTrend = T,
		per_branch_parameters=T, log_anc_states=T)




# plot results
plot_results(fbm_obj, logfile = output_file)
plot_time_varying_trend(fbm_obj, output_file, resfile="trends.pdf")

library(corHMM)
