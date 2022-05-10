# SETUP MAMMALS
# load FBM library
library(ramify)
source("/Users/dsilvestro/Software/fossilBM/fossilBM_lib.R")
setwd("/Users/dsilvestro/switchdrive/Bruna's PhD files/Mammals/FossilBM/")

t_file= "tree1.tre"
tree = read.tree(t_file)

bm_file = "data_bodysize.csv"
bm_data_tbl = read.csv(bm_file) 
bm_data = as.data.frame(bm_data_tbl[,2])
row.names(bm_data) = bm_data_tbl[,1]

d_file = "states_rec_tree1.csv"
anc_trait = read.csv(d_file)

tip_file = "tip_states_tree1.csv"
tip_trait = read.csv(tip_file)
tip_trait[,1] = 1:dim(tip_trait)[1]
tip_states = argmax(tip_trait[,2:5], rows = TRUE)
tip_state_tbl = cbind(tip_trait[,1], as.numeric(tip_states))

# sample from anc state probs
# tip_trait[,1] = 1:dim(tip_trait)[1]
# anc_trait[,1] = 1:dim(anc_trait)[1] + dim(tip_trait)[1]
# comb_trait = rbind(tip_trait, anc_trait)
# state = argmax(comb_trait[,2:5], rows = TRUE)
# state_tbl = cbind(comb_trait[,1], state)

anc_state_file = "lik_anc_states.csv"
anc_state_trait = read.csv(anc_state_file)
anc_trait[,1] = 1:dim(anc_trait)[1] + dim(tip_trait)[1]
colnames(tip_state_tbl) = colnames(anc_state_trait)

state_tbl = rbind(anc_state_trait, tip_state_tbl)




dist_from_the_root = readRDS("tree1_ dist_from_the_root.rds")

# create fossil bm object
fbm_obj <- read_and_transform_data(tree_obj=tree, 
                                   data_obj=bm_data,
                                   log_trait_data=10, 
                                   dist_from_the_root=dist_from_the_root,
                                   state_tbl=state_tbl,
                                   partition_file="order_taxa_no_unique.txt", 
                                   drop_na=TRUE)



# saveRDS(fbm_obj, "/Users/dsilvestro/switchdrive/Bruna's PhD files/Mammals/FossilBM/fbm_obj.rds")
#
#
#
# # start analysis
# source("/Users/dsilvestro/Software/fossilBM/fossilBM_lib.R")
# fbm_obj = readRDS("/Users/dsilvestro/switchdrive/Bruna's PhD files/Mammals/FossilBM/fbm_obj.rds")

fbm_obj$data <- fbm_obj$data / sd(fbm_obj$data, na.rm=T)
fbm_obj$data <- fbm_obj$data - mean(fbm_obj$data, na.rm=T) 
# hist(fbm_obj$data)

output_file <- "/Users/dsilvestro/switchdrive/Bruna's PhD files/Mammals/FossilBM/fbm_mcmc_test_rescalemu00.log"
fbm_obj <- run_mcmc(fbm_obj, logfile=output_file, ngen = 10000, sample = 10, print_freq=10, 
		useTrend = T, linTrend = F,
		per_branch_parameters=F, log_anc_states=F,
        update_mu0=c(), # e.g. c(2, 3, 4) to keep first one const
        estimate_HP=TRUE)






# try plot
fbm_obj = readRDS("/Users/dsilvestro/switchdrive/Bruna's PhD files/Mammals/FossilBM/fbm_obj.rds")



output_file <- "/Users/dsilvestro/switchdrive/Bruna's PhD files/Mammals/FossilBM/fbm_mcmc_test_plot.log"
fbm_obj <- run_mcmc(fbm_obj, logfile=output_file, ngen = 10, sample = 1, print_freq=10, 
		useTrend = T, linTrend = F,
		per_branch_parameters=T, log_anc_states=T)


plot_results(fbm_obj, logfile = output_file)




### OTHER STUFF
# define input files
setwd("/Users/dsilvestro/Software/fossilBM/dev") # set working directory
treefile <- "pruned_tree.tre"
tindex <- 1 # specify which tree should be analyzed if the tree file includes multiple trees
datafile <- "platyrrhine_bm.txt"

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


