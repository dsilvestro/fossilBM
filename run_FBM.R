# define input files
setwd("path") # set working directory
treefile <- "platyrrhine_FBD.trees"
tindex <- 1 # specify which tree should be analyzed if the tree file includes multiple trees
datafile <- "platyrrhine_bodymass.txt"

# load FBM library
source("fossilBM_lib.R")

# process and match input data 
fbm_obj <- read_and_transform_data(treefile, datafile, tindex=tindex, log_trait_data=10)
# Optional arguments:
# log_trait_data = 0 # if set to 0 trait is not log trnsformed
# drop_na = FALSE # if set to TRUE NAs will be dropped otherwise they are inputed
# rescale_trait_data = 1 # multiplier to rescale the trait data
# root_calibration = c(0,100) # normal prior (mean,sd) on root state



# to add a partition file you can use the following:
# fbm_obj <- read_and_transform_data(treefile, datafile, log_trait_data=10, partition_file="clade_partitions.txt")
# this will assume independent rate and trends parameters for each partition (rates and trends are costant within)

# run MCMC
output_file <- "fbm_mcmc.log"
run_mcmc(fbm_obj, logfile=output_file, ngen = 100000, sample = 2500, print_freq=1000)

# Optional arguments:
# useTrend = T # set to FALSE to run a BM model with no trend
# constRate = F # set to TRUE to run a BM model with constant 
# linTrend = F # set to TRUE to run a BM model in which the trend
               # paramater itself follows a linear trend



# plot results
plot_results(fbm_obj, logfile = output_file)
plot_time_varying_trend(fbm_obj, output_file, resfile="trends.pdf")


