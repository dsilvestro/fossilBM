library(TreeSim)

# define input files
setwd("/Users/dsilvestro/Software/fossilBM/") # set working directory

# load FBM library
source("fossilBM_lib.R")

# simulate tree
res_summary = NULL
n_sim = 100

col_names = c("sim_n","ntips", "sigma2", "mu0", "a0", 
			  "BM_sig","BM_sig_m","BM_sig_M","BM_mu0","BM_mu0_m","BM_mu0_M",
			  "BM_a0","BM_a0_m","BM_a0_M","BM_root","BM_root_m","BM_root_M","BM_anc_r2","BM_anc_mse",
		  
			  "BMT_sig","BMT_sig_m","BMT_sig_M","BMT_mu0","BMT_mu0_m","BMT_mu0_M",
			  "BMT_a0","BMT_a0_m","BMT_a0_M","BMT_root","BMT_root_m","BMT_root_M","BMT_anc_r2","BMT_anc_mse",
		  
			  "BMVT_sig","BMVT_sig_m","BMVT_sig_M","BMVT_mu0","BMVT_mu0_m","BMVT_mu0_M",
			  "BMVT_a0","BMVT_a0_m","BMVT_a0_M","BMVT_root","BMVT_root_m","BMVT_root_M","BMVT_anc_r2","BMVT_anc_mse"
			  )

cat(c(col_names, "\n"), sep="\t",file="BMVT_sim.log")

for (sim_i in 1:n_sim){
	
	ntips    = sample(50:150, size=1)
	sigma2   = runif(1, 0.01, 0.1)   # 0.1
	mu0      = runif(1, 0.2, 0.5) * sample(c(-1,1),size=1)   # -1
	a0       = runif(1, 0.02, 0.04) * sample(c(-1,1),size=1) # 0.03

	tree <- sim.bd.taxa(n=ntips, numbsim=1, lambda=0.4*1, mu=0.2*1, frac = 1)[[1]]
	# drop random extinct species
	drop_extinct(tree, rho=0.2)


	# init trait dataframe 
	data = matrix(1, nrow = length(tree$tip.label), ncol = 1)
	rownames(data) = tree$tip.label
	# init 
	fbm_obj <- read_and_transform_data(tree_obj=tree, data_obj=data, log_trait_data=0)


	sim_states <- simulate_trait_data(fbm_obj,  sigma2=sigma2, mu0=mu0, a0=a0, plot=T)
	data <- sim_states[1:fbm_obj$ntips]
	# re-order tip values
	names(data) <- fbm_obj$tree$tip.label 
	data <- data[names(fbm_obj$data)]
	fbm_obj$data <- data
	mu_time_0 = mu0 + a0*max(fbm_obj$dist_from_midpoint)
	mu_root = mu0 + a0*min(fbm_obj$dist_from_midpoint)


	# run MCMCs
	model = "BM"
	output_file <- paste0("fbm_mcmc", model)
	run_mcmc(fbm_obj, logfile=paste0(output_file, ".log"), 
			 useTrend = F, constRate = T, linTrend = F, bdmcmc_freq = 0,
			 ngen = 2000, sample = 40, print_freq=1000)

	x0 = plot_results(fbm_obj, logfile = paste0(output_file, ".log"), resfile=paste0(output_file, ".pdf"),return_estimated_prm=T)
	# mean(abs(x0-sim_states)/mean(sim_states))


	model = "BMT"
	output_file <- paste0("fbm_mcmc", model)
	run_mcmc(fbm_obj, logfile=paste0(output_file, ".log"), 
			 useTrend = T, constRate = T, linTrend = F, bdmcmc_freq = 0,
			 ngen = 2000, sample = 40, print_freq=1000)

	x1 = plot_results(fbm_obj, logfile = paste0(output_file, ".log"), resfile=paste0(output_file, ".pdf"),return_estimated_prm=T)

	model = "BMVT"
	output_file <- paste0("fbm_mcmc", model)
	run_mcmc(fbm_obj, logfile=paste0(output_file, ".log"), 
			 useTrend = T, constRate = T, linTrend = T, bdmcmc_freq = 0,
			 ngen = 5000, sample = 100, print_freq=1000)

	x2 = plot_results(fbm_obj, logfile = paste0(output_file, ".log"), resfile=paste0(output_file, ".pdf"),return_estimated_prm=T)

	if (2<1){
		par(mfrow=c(1,3))
		plot(sim_states, x0[[1]])
		abline(coef=c(0,1), lty=2)

		plot(sim_states, x1[[1]])
		abline(coef=c(0,1), lty=2)

		plot(sim_states, x2[[1]])
		abline(coef=c(0,1), lty=2)
	}



	r2_mse_BM <- get_r2_mse(sim_states, x0[[1]])
	r2_mse_BMT <- get_r2_mse(sim_states, x1[[1]])
	r2_mse_BMVT <- get_r2_mse(sim_states, x2[[1]])


	res_summary_row = c(sim_i, ntips, sigma2, mu0, a0, 
						x0[[2]], r2_mse_BM, x1[[2]], 
						r2_mse_BMT, x2[[2]], r2_mse_BMVT )
						
	names(res_summary_row) = col_names
	res_summary <- rbind(res_summary, res_summary_row)
	cat(c(res_summary_row, "\n"), sep="\t",file="BMVT_sim.log", append=T)
	







	# res_all_files = sprintf("%s.rda", filename)
	# save.image(res_all_files)

	
}


res_summary_df = as.data.frame(res_summary,row.names=F)
colnames(res_summary_df) = col_names

