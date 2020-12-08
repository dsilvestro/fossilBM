# define input files
setwd("wd") # set working directory
treefile <- "example_files/platyrrhine_FBD.trees"
tindex <- 1 # specify which tree should be analyzed if the tree file includes multiple trees
datafile <- "example_files/platyrrhine_bodymass.txt"

# load FBM library
source("fossilBM_lib.R")

# process and match input data 
fbm_obj <- read_and_transform_data(treefile, datafile, log_trait_data=0, drop_na=F)



simulate_trait_data <- function(fbm_obj, sigma2=0.2, mu0=0, a0=0){
	ntips	  <- fbm_obj$ntips
	tree <- fbm_obj$tree
	D		  <- fbm_obj$D
	prior_tbl  <- fbm_obj$prior_tbl
	root_dist <- fbm_obj$dist_from_the_root
    dist_from_midpoint <- fbm_obj$dist_from_the_root
	
	root_state = 0
	all_states = rep(root_state, (ntips*2-1))
	
	for (indx in (ntips+1):(ntips*2-1)){
		
		i = which(D[,1]==indx)
		
		anc_ind <- D[i,1];
		a_ind   <- D[i,2];  # index of descendants
		b_ind   <- D[i,3];  # index of descendants
		vpa	 <- D[i,4];  # br length
		vpb	 <- D[i,5];  # br length
		anc = all_states[anc_ind]

		m1 <- anc + (vpa*mu0 + a0*vpa*dist_from_midpoint[a_ind])
		s1 <- sqrt((vpa)*sigma2)
		m2 <- anc + (vpb*mu0 + a0*vpb*dist_from_midpoint[b_ind])
		s2 <- sqrt((vpb)*sigma2)
		
		all_states[a_ind] <- rnorm(1, m1, s1)
		all_states[b_ind] <- rnorm(1, m2, s2)
		
	}
	plot(-(max(root_dist)-root_dist), all_states, pch=16)
    points(-(max(root_dist)-root_dist)[1:fbm_obj$ntips], all_states[1:fbm_obj$ntips], pch=16, col="red")
	return(all_states)
}	
	

sigma2=0.1
mu0= -1
a0=0.03
sim_states <- simulate_trait_data(fbm_obj,  sigma2=sigma2, mu0=mu0, a0=a0)
data <- sim_states[1:fbm_obj$ntips]
# re-order tip values
names(data) <- fbm_obj$tree$tip.label 
data <- data[names(fbm_obj$data)]
fbm_obj$data <- data
mu_time_0 = mu0 + a0*max(fbm_obj$dist_from_midpoint)
mu_root = mu0 + a0*min(fbm_obj$dist_from_midpoint)


a0	8.884E-3	<html><font color="#EE0000">13</font></html> 	R



# run MCMCs
model = "BM"
output_file <- paste0("fbm_mcmc", model)
run_mcmc(fbm_obj, logfile=paste0(output_file, ".log"), 
		 useTrend = F, constRate = T, linTrend = F, bdmcmc_freq = 0,
		 ngen = 15000, sample = 50, print_freq=1000)

x0 = plot_results(fbm_obj, logfile = paste0(output_file, ".log"), resfile=paste0(output_file, ".pdf"),return_anc_states=T)
# mean(abs(x0-sim_states)/mean(sim_states))


model = "BMT"
output_file <- paste0("fbm_mcmc", model)
run_mcmc(fbm_obj, logfile=paste0(output_file, ".log"), 
		 useTrend = T, constRate = T, linTrend = F, bdmcmc_freq = 0,
		 ngen = 15000, sample = 50, print_freq=1000)

x1 = plot_results(fbm_obj, logfile = paste0(output_file, ".log"), resfile=paste0(output_file, ".pdf"),return_anc_states=T)

model = "BMVT"
output_file <- paste0("fbm_mcmc", model)
run_mcmc(fbm_obj, logfile=paste0(output_file, ".log"), 
		 useTrend = T, constRate = T, linTrend = T, bdmcmc_freq = 0,
		 ngen = 50000, sample = 50, print_freq=1000)

x2 = plot_results(fbm_obj, logfile = paste0(output_file, ".log"), resfile=paste0(output_file, ".pdf"),return_anc_states=T)

par(mfrow=c(1,3))
plot(sim_states, x0)
abline(coef=c(0,1), lty=2)

plot(sim_states, x1)
abline(coef=c(0,1), lty=2)

plot(sim_states, x2)
abline(coef=c(0,1), lty=2)

