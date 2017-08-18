#!/usr/bin/Rscript
# version 20150713


# module add R/3.0.2
# cd /scratch/ul/monthly/dsilvest/BM && RScript mcmcFossilBM-BD.R --wd /scratch/ul/monthly/dsilvest/BM --j %s --v 8 --s 2 --r 0.1 --t 100


arg <- commandArgs(trailingOnly=TRUE)
library(optparse)
library(scales)

use_loc_libraries = 0
useBMT=TRUE 

require(phytools)
library(geiger)
library(TreeSim)

##################################################

option_list <- list(                   
	make_option("--n",  type="integer",         default=250000, help=("MCMC iterations (%default)"),metavar="iterations"),
	make_option("--wd", type="character",       default="",     help=("working dir"),metavar="wd"),
	make_option("--plot_res", type="integer",   default=1,    help=("plot pdf file (%default)"), metavar="plotf"),

	# SIMULATION SETTINGS
	make_option("--t",    type="integer",  default=50, help=("Tree size (%default)"),metavar="ntips"),
	make_option("--v",    type="double",   default=16, help=("Magnitude rate shift (%default)"),metavar="xfold"),    
	make_option("--r",    type="double",   default=0,  help=("Baseline rate (%default)"),metavar="rate"),    
	make_option("--s",    type="integer",  default=0,  help=("no. shifts sig2 (%default)"),metavar="shifts"),
	make_option("--mu",   type="double",   default=NA, help=("mu0 (%default)"),metavar="shifts"),
	make_option("--sm",   type="integer",  default=0,  help=("no. shifts mu0 (%default)"),metavar="shifts"),
	make_option("--rda",  type="integer",  default=0,  help=("if 1 use existing RDA file (%default)"), metavar="rda"),
	make_option("--qfos", type="integer",  default=20, help=("Mean number of simulated fossils"), metavar="qfos"),
	make_option("--j",    type="integer",  default=1,  help=("Simulation number (%default)"),metavar="replicate"),
	
	# COMMANDS FOR EMPIRICAL ANALYSIS
	make_option("--out",      type="character", default="",   help=("name of (%default)"), metavar=""),
	make_option("--tfile",    type="character", default="NA", help=("Name of tree nexus file"), metavar="tfile"),
	make_option("--tindex",   type="integer",   default=1,    help=("Index of the tree"), metavar="tindex"),
	make_option("--dfile",    type="character", default="NA", help=("Name of trait file"), metavar="dfile"),
	make_option("--rmF",      type="integer",   default=0,    help=("Set to 1 to remove extinct tips"), metavar="rmF"),
	make_option("--log",      type="double",    default=0,    help=("log transform data - opts: 10 (log10), 1 (ln), 0 (no transform)"),metavar="log"),
	make_option("--rescale",  type="double",    default=1,    help=("multiply trait"),metavar="log")
	)

parser_object <- OptionParser(usage = "Usage: %prog [Options]", option_list=option_list, description="...")
opt <- parse_args(parser_object, args = commandArgs(trailingOnly = TRUE), positional_arguments=TRUE)
	
ntips              = opt$options$t  
nGibbs_gen         = opt$options$n
Gibbs_sample       = 250
root_calibration   = c(0,100)
print_f            = 250
xfold              = opt$options$v
dynamicPlot        = F
w_dir              = opt$options$wd 
rep_number         = opt$options$j
sim_n              = 1
sig2               = opt$options$r
n_shifts           = opt$options$s
mu0                = opt$options$mu
n_shifts_mu0       = opt$options$sm
plot_res           = opt$options$plot_res
output	       = opt$options$out
q_f		       = opt$options$qfos
treefile	       = opt$options$tfile
tindex	       = opt$options$tindex
dfile		       = opt$options$dfile
remF 		       = opt$options$rmF
log_trait_data     = opt$options$log
rescale_trait_data = opt$options$rescale

# qfos = 20 is the default, qfos can be any number e.g. 5, but if qfos is 1, it takes only the oldest fossil of ~20 simulated

if (sig2==0){
	sig2 = round(rgamma(1,2,5),2)
}

if (is.na(mu0)){
	mu0 = rnorm(1, 2,0.5)*sample(c(-1,1),1)
}

setwd(w_dir)

use_RDA = opt$options$rda
extant_only=0 # set to 1 to run on exant only

if (use_RDA==1){
	res_all_files = sprintf("sim_%s%s.rda", opt$options$j, opt$options$out)
	load(res_all_files)
	rda_exist=1
	print(paste("Loading ", res_all_files))
	use_RDA=1
	print(use_RDA)
}


update_multiplier_proposal <- function(i,d){
	u = runif(1,0,1)
	l = 2*log(d)
	m = exp(l*(u-0.5))
 	ii = i * m
	U=log(m)
	#if (ii>10){
	#	return(c(i,0))
	#}else{return(c(ii, U))}
	return(c(ii, U))
}

fN<-function(xa, xb, vpa, vpb, sigma2) {
	#the normal density. Sigma2 is the variance. Same as log(dnorm(xa-xb, 0, sqrt(sigma * sigma * (vpa + vpb))))
	return( dnorm((xa-xb),mean=0,sd= sqrt(sigma2 * (vpa + vpb)), log=T) )
}

fN2<-function(xa, xb, vpa, vpb, sigma2_a,sigma2_b, anc) {
	#same as fN but with a known ancestor instead of xa-xb
	return(dnorm(xa, mean=anc, sd= sqrt(vpa*sigma2_a),log=T) + dnorm(xb, mean=anc, sd= sqrt(vpb*sigma2_b),log=T)  )
}

newlnLike<-function(tree, vector_tip_root_nodes_values, sigma2,D,mu0) {
	vec_values = vector_tip_root_nodes_values
	anc_ind_v <- D[,1]
	a_ind_v   <- D[,2]
	b_ind_v   <- D[,3]
	vpa_v     <- D[,4]
	vpb_v     <- D[,5]
	anc_v <- vec_values[anc_ind_v]
	a_v   <- vec_values[a_ind_v]
	b_v   <- vec_values[b_ind_v]
	a_ind_v[a_ind_v>ntips] =a_ind_v[a_ind_v>ntips]-1
	b_ind_v[b_ind_v>ntips] =b_ind_v[b_ind_v>ntips]-1
	L1 = dnorm(a_v, mean=(anc_v+mu0[a_ind_v]*vpa_v), sd= sqrt(vpa_v*sigma2[a_ind_v]),log=T)
	L2 = dnorm(b_v, mean=(anc_v+mu0[b_ind_v]*vpb_v), sd= sqrt(vpb_v*sigma2[b_ind_v]),log=T)  
	L = L1+L2
	return(L) 
}

build_table <- function (tree,ntips,data){
	tree.edge.ordered<-tree$edge[order(tree$edge[,2]),]
	# the root is at ntips+1
	root<-ntips+1
	# the table with ancestor, the two descendants and the two branch length
	table.tree<-matrix(NA,ntips-1,5)
	colnames(table.tree)<-c("ancestor","descendant1","descendant2","branch.length1", "branch.length2")

	#contains the root and all the nodes 
	table.tree[,1]<-c(ntips+1,tree.edge.ordered[(ntips+1):(ntips+ntips-2),2])

	for(i in 1:dim(table.tree)[1]){
		table.tree[i,c("descendant1","descendant2")] <- tree$edge[which(tree$edge[,1]==table.tree[i,1]),2]
		table.tree[i,c("branch.length1","branch.length2")] <- tree$edge.length[which(tree$edge[,1]==table.tree[i,1])]
	}
	return(table.tree)
}

get_calibration_tbl <- function (D,root_calibration){
	calibration_tbl =NULL	
	for (i in 1:dim(D)[1]){
		calibration_tbl = rbind(calibration_tbl,c(0,100))
	}	
	row.names(calibration_tbl)=D[,1]
	calibration_tbl[1,] = root_calibration
	return(calibration_tbl)
}

################ RUN GIBBS FUNCTIONS
get_joint_mu_s2 <- function (mu_f,s2_f,mu_g,s2_g){
	s2_fg = (s2_f*s2_g)/(s2_f+s2_g)
	mu_fg = (mu_f/s2_f + mu_g/s2_g) * s2_fg
	return(c(mu_fg,s2_fg))
}

runGibbs <- function(sigma2, vector_tip_root_nodes_values,D,prior_tbl,mu0,get_expected=0) {
	#prior_tbl = calc_prior_dist_root_calibration(D,sigma2)
	
	vec_values = vector_tip_root_nodes_values
        # loop over Johnatan's tbl from most recent to root
        for (i in dim(D)[1]:2){
      	anc_ind <- D[i,1];
      	a_ind   <- D[i,2];  # index of descendants
      	b_ind   <- D[i,3];  # index of descendants
      	vpa     <- D[i,4];  # br length
      	vpb     <- D[i,5];  # br length
      	anc = vec_values[anc_ind]
      	a = vec_values[a_ind]
      	b = vec_values[b_ind]	
		# calibration prior
		calibrated_prior_mu = prior_tbl[which(rownames(prior_tbl)==anc_ind),1]
		calibrated_prior_s2 = prior_tbl[which(rownames(prior_tbl)==anc_ind),2]
		
		if (a_ind>ntips){
			sig_a_ind=a_ind-1 # because no sig2 value for root
		}else{sig_a_ind=a_ind}

		if (b_ind>ntips){
			sig_b_ind=b_ind-1
		}else{sig_b_ind=b_ind}

		desc_lik = get_joint_mu_s2((a-vpa*mu0[sig_a_ind]), (vpa)*sigma2[sig_a_ind], (b-vpb*mu0[sig_b_ind]), (vpb)*sigma2[sig_b_ind])
		prior_prm = prior_tbl[which(rownames(prior_tbl)==anc_ind),]
		
		# lik from (stem) ancestral state
		if (i>1){
			if (anc_ind %in% D[,2]){
				stem_ind = D[which(D[,2]==anc_ind),1]
				stem_brl = D[which(D[,2]==anc_ind),4]
			}
			else{
				stem_ind = D[which(D[,3]==anc_ind),1]
				stem_brl = D[which(D[,3]==anc_ind),5]	
			}
			
			sig_stem_ind=stem_ind-1
			stem_val = vec_values[stem_ind]		
			stem_lik = c((stem_val + stem_brl*mu0[sig_stem_ind]) , (stem_brl)*sigma2[sig_stem_ind])
			
			lik_val = get_joint_mu_s2(stem_lik[1],stem_lik[2],desc_lik[1],desc_lik[2])

			post_prm= get_joint_mu_s2(lik_val[1],lik_val[2],prior_prm[1], prior_prm[2])
			if (get_expected==1){
				vec_values[anc_ind] = post_prm[1]
			}else{
				vec_values[anc_ind] = rnorm(1, mean=post_prm[1], sd=sqrt(post_prm[2]))
			}
		}
	}
	return(vec_values[D[,1]] )
}

calc_prior <- function(sig2, a, y, mu0) {
	prior_sig2 = sum(dexp(sig2, 0.5, log = TRUE) ) #sum(dgamma(sig2, 2,1, log = TRUE) )
	prior_root = sum(dnorm(c(a), mean = prior_tbl[,1], sd = prior_tbl[,2], log = T))
	prior_anc  = sum(dnorm(c(a, y), mean = 0, sd = 100, log = T))
	prior_mu0  = sum(dnorm(mu0, mean = 0, sd = 1, log = T))
	return(prior_sig2+prior_root+prior_anc+prior_mu0)
}

################ BDMCMC SAMPLER

calc_rel_prob <- function(log_lik){
	rel_prob=exp(log_lik-max(log_lik))
	return (rel_prob/sum(rel_prob)	)
}

random_choice_P <-function(vector){
	probDeath=cumsum(vector/sum(vector)) # cumulative prob (used to randomly sample one 
	r=runif(1,0,1)                          # parameter based on its deathRate)
	probDeath=sort(append(probDeath, r))
	ind=which(probDeath==r)
	return (c(vector[ind], ind))
}

calc_deathRates <- function (tree, trait_states, sig2, ind_sig, mu0, ind_mu0,D,IND_edge,desc_index_list,death_sig2=1){
	list_deathrates=c(0) # can't kill first
	current_lik = sum(newlnLike(tree, trait_states, sig2[ind_sig],D,mu0[ind_mu0]))
	if (death_sig2==1){
		N=length(unique(ind_sig))
	}else{N=length(unique(ind_mu0))}
	
	for (i in 2:N){
		if (death_sig2==1){
			multi_ind = ind_sig
		}else{
			multi_ind = ind_mu0
		}
		multi_ind_edge  = multi_ind[IND_edge]
		# get edges with index i
		Yi = which(multi_ind_edge==i)
		crown_node_ind_i = min(tree$edge[Yi,1])
		if (crown_node_ind_i==(ntips+1)){ # ROOT
			ind_of_ancestor=1
		} else{
			ind_of_ancestor= multi_ind_edge[which(tree$edge[,2]==crown_node_ind_i)]
		}

		if (death_sig2==1){
			new_ind_sig = ind_sig
			new_ind_sig[new_ind_sig==i] = ind_of_ancestor		
			rates_temp= sig2[new_ind_sig]
			mu0_temp = mu0[ind_mu0]
		}else{
			new_ind_mu0 = ind_mu0
			new_ind_mu0[new_ind_mu0==i] = ind_of_ancestor		
			rates_temp= sig2[ind_sig]
			mu0_temp = mu0[new_ind_mu0]
		}
		
		list_deathrates[i] = sum(newlnLike(tree, trait_states, rates_temp,D,mu0_temp))
	}
	list_deathrates = exp(list_deathrates - current_lik)	
	list_deathrates[1] = 0
	return(list_deathrates)
}

run_BDMCMC <- function(tree, trait_states, sig2, ind_sig, mu0, ind_mu0,D,IND_edge,new_ind,desc_index_list,death_sig2=1){
	init_lik = sum(newlnLike(tree, trait_states, sig2[ind_sig],D,mu0[ind_mu0]))
	cont_time=0
	len_cont_time=1
	birthRate=1
	# cat("\n\nSTART\n")
	# print(c("initial rates",sig2))
	while (cont_time<len_cont_time){
		if (death_sig2==1){
			deathRate = rep(0, length(unique(ind_sig)))
			n_likBD   = rep(0, length(unique(ind_sig)))
		}else{
			deathRate = rep(0, length(unique(ind_mu0)))
			n_likBD   = rep(0, length(unique(ind_mu0)))
		}
		
		if (length(unique(ind_sig))>1 && death_sig2==1){
				deathRate = calc_deathRates(tree, trait_states, sig2, ind_sig, mu0, ind_mu0, D,IND_edge,desc_index_list,1)
		}else if (length(unique(ind_mu0))>1 && death_sig2==2){
				deathRate = calc_deathRates(tree, trait_states, sig2, ind_sig, mu0, ind_mu0, D,IND_edge,desc_index_list,2)
		}
		
		deltaRate=sum(deathRate)
		
		cont_time = cont_time + rexp(1, rate = min((deltaRate+birthRate), 100000))
		
		if ((deltaRate+birthRate)>100 && cont_time>len_cont_time){
			cont_time <- len_cont_time
		}
		
		if (cont_time>len_cont_time){
			"pass"			
		}else{
			Pr_birth= birthRate/(birthRate+deltaRate)
			Pr_death= 1-Pr_birth
			
			if (runif(1)<Pr_birth){ # ADD PARAMETER
				rates_temp= sig2[ind_sig]
				trend_temp= mu0[ind_mu0]
				
				if (death_sig2==1){
					n_ind  = sample(1:length(ind_sig),1)
					multi_ind = ind_sig
				}else{
					n_ind  = sample(1:length(ind_mu0),1)
					multi_ind = ind_mu0					
				}
				
				multi_ind_edge  = multi_ind[IND_edge]
				multi_ind_new = multi_ind_edge
				# index that is being affected
				split_ind = multi_ind_new[n_ind]
				
				ind_desc = desc_index_list[[n_ind]]
				if (length(ind_desc) ==1){ # NO shifts allowed on single branches (tips)
					"pass"
				}else{
				
					if (death_sig2==1){
						for (j in ind_desc){
							if (multi_ind_new[j]==split_ind){
								multi_ind_new[j] = max(ind_sig)+1
							}
						}
						# BDMCMC PROPOSAL
						sig2_original  =  sig2
						ind_sig_original  = ind_sig
						sig2=c(sig2,rexp(1,0.5)) # sample from rate prior
						# ind_sig in the original order
						ind_sig=multi_ind_new[new_ind]
						if (length(sig2)>length(unique(ind_sig))){
							# the new index overlaps completely with the previous
							# remove unused sig2
						 	sig2= sig2_original
							ind_sig= ind_sig_original
						}
					}else{
						for (j in ind_desc){
							if (multi_ind_new[j]==split_ind){multi_ind_new[j] = max(ind_mu0)+1}
						}
						# BDMCMC PROPOSAL
						mu0_original  =  mu0
						ind_mu0_original  = ind_mu0
						mu0=c(mu0,rnorm(1,0,1)) # sample from rate prior
						# ind_sig in the original order
						ind_mu0=multi_ind_new[new_ind]
						if (length(mu0)>length(unique(ind_mu0))){
							# the new index overlaps completely with the previous
							# remove unused mu0
						 	mu0= mu0_original
							ind_mu0= ind_mu0_original
						}
					} 
					
				}
			}else{ # REMOVE PARAMETER
					i=random_choice_P(deathRate)[2]
				
					if (death_sig2==1){
						multi_ind = ind_sig
					}else{
						multi_ind = ind_mu0
					}
					
					
					multi_ind_edge  = multi_ind[IND_edge]
					# get edges with index i
					Y = which(multi_ind_edge==i)
					#__ cat("\ndead edge:",Y, i)
					crown_node_ind_i = min(tree$edge[Y,1])
				
					if (crown_node_ind_i<Inf && crown_node_ind_i>-Inf){
						"ok"
					}else{
						cat("\ndead crown:",crown_node_ind_i, i,"\n")
						cat("\ndead edge:",Y,"ind:" ,i,"\nY:",tree$edge[Y,1])
					}
				
					if (crown_node_ind_i==(ntips+1)){ # ROOT
						ind_sig_of_ancestor_i=1
					}else{
						ind_sig_of_ancestor_i= multi_ind_edge[which(tree$edge[,2]==crown_node_ind_i)]
					}
				
					if (death_sig2==1){
						# update indexes
						new_ind_sig = ind_sig
						new_ind_sig[new_ind_sig==i] = ind_sig_of_ancestor_i
						ind_sig = new_ind_sig		
						# remove unused sig2
						sig2 = sig2[-i]				
						# rescale indexes
						ind_sig[ind_sig>i]=ind_sig[ind_sig>i]-1
					}else{
						# update indexes
						new_ind_mu0 = ind_mu0
						new_ind_mu0[new_ind_mu0==i] = ind_sig_of_ancestor_i
						ind_mu0 = new_ind_mu0		
						# remove unused sig2
						mu0 = mu0[-i]				
						# rescale indexes
						ind_mu0[ind_mu0>i]=ind_mu0[ind_mu0>i]-1					
					}

			}
			
		}
		
	}
	#__ print(list(sig2,ind_sig))	
	return (list(sig2,ind_sig,mu0,ind_mu0))
}

################################## START MCMC ###################################
mcmc.gibbs4 <- function (tree, x, D, prior_tbl,true_rate,ngen = 100000, control = list(),useVCV=F, sample=100,logfile="log",update_sig_freq=0.5,dynamicPlot = F,useFixPart=F,bdmcmc_freq=0.75,useTrend=F,print_freq=100){
	x=data
	TE=as.matrix(tree$edge)
	cat(c("it", "posterior","likelihood","prior", "sig2", "mu0", "MAPE","sdAPE","K_sig2","K_mu0", paste("sig2", 1:length(TE[,1]),sep="_"),paste("mu0", 1:length(TE[,1]),sep="_"), paste("anc",D[,1],sep="_"),"\n"),sep="\t", file=logfile, append=F)
	a <- 0
	y <- rep(0, tree$Nnode - 1)
	sig2 <- c(0.2)  
	mu0 <- 0
	ind_sig2 <- rep(1, length(c(x, y)))
	ind_mu0  <- rep(1, length(c(x, y)))
	mean_sig <- rep(0, length(c(x, y)))
	mean_anc <- rep(0, length(c(a, y)))
	if (dynamicPlot == T){dev.new()}
	
	#names(ind_sig2) = c(names(x), D[-1,1])

	x <- x[tree$tip.label]
	if (is.null(names(y))) {
		names(y) <- length(tree$tip) + 2:tree$Nnode
	}
	else {y[as.character(length(tree$tip) + 2:tree$Nnode)]}

	L  <- newlnLike(tree, c(x, a, y), sig2[ind_sig2],D,mu0[ind_mu0])
	Pr <- calc_prior(sig2, a, y,mu0)


	# get indexes
	IND_edge = c()
	for (ind_edge in tree$edge[,2]){
		if (ind_edge>ntips){
			ind_edge_t=ind_edge-1
		}else{ind_edge_t=ind_edge}
		IND_edge = append(IND_edge,ind_edge_t)
	}
	# indexes sorted by edge index
	multi_sig = rep(0,length(IND_edge))
	names(multi_sig) = 1:length(multi_sig)
	multi_sig_edge  = multi_sig[IND_edge]
	new_ind=order(as.numeric(names(multi_sig_edge)))

	# list index descendants
	desc_index_list = list()
	for (i in 1:length(multi_sig)){
		node_ind=tree$edge[i,2]
		ind_desc_nodes = getDescendants(tree, node=node_ind)
		ind_desc_edges = which(tree$edge[,2]%in% ind_desc_nodes == T)
		V = unique(c(i,ind_desc_edges))
		desc_index_list[[i]] = unique(c(i,ind_desc_edges))
	}

	# START MCMC
	for (i in 1:ngen) {
    
		y.prime     = y
		L.prime     = L
		Pr.prime    = Pr
		a.prime     = a
		sig2.prime  = sig2
		mu0.prime   = mu0
		gibbs=0
		hastings=0
    	    	
		# SCREEN OUTPUT / DYNAMIC PLOTS
		if (i %% print_freq==0){
			cat(c("\n",i,round(sum(L),2)))
		}
		
		if (dynamicPlot == T){
			mean_sig = rbind(mean_sig,sig2[ind_sig2])
			#mean_anc = rbind(mean_anc,c(x,a,y))
			mean_anc = rbind(mean_anc,c(a,y))
		
			if (i%%print_f == 0 || i==1) {
				cat(c("\n",i,round(c(sum(L),sig2,mu0),2)))
				par(mfrow=c(1,2))
				start = round(dim(mean_sig)[1]*0.1)
				end = dim(mean_sig)[1]
	    		        #rates_temp =sig2.prime[ind_sig2]
	    		        rates_temp = apply(mean_sig[start:end,], 2, FUN=median)
	    		        plot.phylo(tree, edge.width=rates_temp[IND_edge]*3, main="median rates",show.tip.label = F)
				
	    		        rates_temp =sig2.prime[ind_sig2]
	    		        plot.phylo(tree, edge.width=rates_temp[IND_edge]*3, main=paste("lik:",round(sum(L),2),"cat:",length(unique(ind_sig2))),show.tip.label = F)
				
				#phenogram(tree, c(x,a,y),ylim=c(-1,1)) #,add=add,col=col	
				# delta anc-states
				#start = round(dim(mean_anc)[1]*0.1)
				#end = dim(mean_anc)[1]
	    		        #anc_temp = apply(mean_anc[start:end,], 2, FUN=mean)
				#edge_w=abs(anc_temp[tree$edge[,1]] - anc_temp[tree$edge[,2]])/sqrt(tree$edge.length)
				### hist error
				##Err = (rates_temp[IND_edge]-true_rate)
				##hist(Err, main = paste("MAE:",round(mean(abs(Err)),3), "Qtl:", round(sort(Err)[as.integer(length(Err)*0.25)],3), round(sort(Err)[as.integer(length(Err)*0.75)],3) )  ) #, "MAPE:",round(mean(abs(anc_temp-true_rate)true_ratec),3)) ) 
				
				
				
			}
		}
    
		j <- (i - 1)%%(tree$Nnode + 1)
		rr= runif(2,0,1)
		
		update_sig_freq=0.5
		if (rr[1]<update_sig_freq) {
			if (rr[2]<0.33){ # SIG2 UPDATE
				s_ind  = sample(1:length(sig2),1)
				sig2_update <-  update_multiplier_proposal(sig2.prime[s_ind],1.2)
				sig2.prime[s_ind] = sig2_update[1]
				hastings = sig2_update[2]			
			}
			else if (rr[2]<0.66 && useTrend==T){ # MU0 UPDATE
				m_ind  = sample(1:length(mu0),1)
				mu0_update <-  mu0[m_ind] + rnorm(n = 1, sd = 0.5)
				mu0.prime[m_ind] = mu0_update
			}
			else{ # ROOT STATE
				a.prime <- a + rnorm(n = 1, sd = 0.5) #sqrt(con$prop[j + 1]))
			}
		}
		 else{
			# ANC STATES
			RUNIF = runif(3,0,1)
			if (RUNIF[1]<bdmcmc_freq && i > 1000){
				trait_states=c(x, a.prime, y.prime)
				if (useTrend==T){
					death_sig2=sample(2,1)
				}else{death_sig2=1}
				bdmcmc_list   = run_BDMCMC(tree, trait_states, sig2.prime,ind_sig2, mu0.prime, ind_mu0,D,IND_edge,new_ind,desc_index_list,death_sig2)
				sig2.prime    = bdmcmc_list[[1]]
				ind_sig2      = bdmcmc_list[[2]]					
				mu0.prime     = bdmcmc_list[[3]]
				ind_mu0       = bdmcmc_list[[4]]										
				#print (ind_mu0)
				gibbs=1
			} else {
			
			vector_tip_root_nodes_values = c(x, a.prime, y.prime)
			y_temp = runGibbs(sig2.prime[ind_sig2], vector_tip_root_nodes_values,D,prior_tbl,mu0.prime[ind_mu0],get_expected=0)
			a.prime= y_temp[1]
			y.prime= y_temp[-1]
			gibbs=1
			}
		}

	       # calc post
		L.prime <- newlnLike(tree, c(x, a.prime, y.prime), sig2.prime[ind_sig2],D,mu0.prime[ind_mu0])	
		Pr.prime <- calc_prior(sig2.prime, a.prime, y.prime, mu0.prime)
	
		if ( (sum(Pr.prime) + sum(L.prime) - sum(Pr) - sum(L) + hastings) >= log(runif(1,0,1)) || gibbs==1){    
			y    = y.prime
			L    = L.prime
			Pr   = Pr.prime
			a    = a.prime
			sig2 = sig2.prime
			mu0  = mu0.prime
		}
 
	     if (i%%sample == 0) {
			rates_temp=sig2[ind_sig2]
			trends_temp=mu0[ind_mu0]
			MAPE = mean(abs((rates_temp[IND_edge]-true_rate))/true_rate)
			sdAPE =  sd(abs((rates_temp[IND_edge]-true_rate))/true_rate)
			cat(c(i,sum(L)+sum(Pr), sum(L),sum(Pr), mean(sig2[ind_sig2]),mean(mu0[ind_mu0]), MAPE, sdAPE, length(sig2),length(mu0), rates_temp[IND_edge],trends_temp[IND_edge], a, y, "\n"),sep="\t", file=logfile, append=T) 
	    }
    
}
}

################################## END   MCMC ###################################

############################### SIMULATE DATA #######################################
simBMT_shifts<-function(tree,a=0,mu=0,sig2=1,bounds=c(-Inf,Inf),internal=T,nsim=1){
	if(bounds[2]<bounds[1]){
		warning("bounds[2] must be > bounds[1]. Simulating without bounds.")
		bounds<-c(-Inf,Inf)
	}
	if(bounds[1]==-Inf&&bounds[2]==Inf) no.bounds=TRUE
	else no.bounds=FALSE
	if(a<bounds[1]||a>bounds[2]){
		warning("a must be bounds[1]<a<bounds[2]. Setting a to midpoint of bounds.")
		a<-bounds[1]+(bounds[2]-bounds[1])/2
	}
	if(length(which(sig2<0))>0){
		warning("sig2 must be > 0.  Setting sig2 to 1.0.")
		sig2[which(sig2<0)]<-1.0
	}
	# function for reflection off bounds
	reflect<-function(yy,bounds){
		while(yy<bounds[1]||yy>bounds[2]){
			if(yy<bounds[1]) yy<-2*bounds[1]-yy
			if(yy>bounds[2]) yy<-2*bounds[2]-yy
		}
		return(yy)
	}
	
	nsp<-length(tree$tip) #nb of species
	
	# first simulate changes along each branch
	nedges<-length(tree$edge.length) #number of edges in the tree
	ndraws<-nedges*nsim #total number of calls to rnorm
	x<-matrix(data=rnorm(n=ndraws,mean=rep(mu*tree$edge.length,nsim),sd=rep(sqrt(sig2*tree$edge.length),nsim)),nedges,nsim)
	# now add them up
	y<-array(0,dim=c(nrow(tree$edge),ncol(tree$edge),nsim))
	for(i in 1:nrow(x)){
		if(tree$edge[i,1]==(nsp+1))
			y[i,1,]<-a
		else
			y[i,1,]<-y[match(tree$edge[i,1],tree$edge[,2]),2,]

		y[i,2,]<-y[i,1,]+x[i,]
		if(!no.bounds) y[i,2,]<-apply(as.matrix(y[i,2,]),1,function(yy) reflect(yy,bounds))
	}
	
	x<-matrix(data=rbind(y[1,1,],as.matrix(y[,2,])), nedges+1,nsim)
	rownames(x)<-c(nsp+1,tree$edge[,2])
	x<-as.matrix(x[as.character(1:(nsp+tree$Nnode)),])
	rownames(x)[1:nsp]<-tree$tip.label
	
	if(internal==TRUE) {
		return(x[1:nrow(x),]) # include internal nodes
	}
	else {
		return(x[1:length(tree$tip.label),]) # tip nodes only
	}
}

sim_data_bmt_shifts <- function(s2,xfold_sig2,n_shifts_sig2,n_shifts_mu0,mu0){

	print("Simulating trees")
	tree<-sim.bd.taxa(n=ntips, numbsim=1, lambda=0.4, mu=0.2, frac = 1)[[1]]
	original_simulated_tree <-tree
	
	### Victor's script for taking the simulated tree, removing extinct sp and adding fossils info according to q 
	
	if (q_f > 1) {
		q_fossil <- q_f/sum(tree$edge.length)
		temp_tree = sprintf("temp_%s_n.tree", rep_number)
		write.tree(tree, temp_tree)
		FBD_temp = sprintf("FBDtemp_%s_n.tree", rep_number)
		system(paste('python simulate_fossils_on_tree.py ',temp_tree, FBD_temp, q_fossil))
		tree <- read.tree(FBD_temp)
		}
		
	if (q_f == 0) {
		q_fossil <- 0
		temp_tree = sprintf("temp_%s_n.tree", rep_number)
		write.tree(tree, temp_tree)
		FBD_temp = sprintf("FBDtemp_%s_n.tree", rep_number)
		system(paste('python simulate_fossils_on_tree.py ',temp_tree, FBD_temp, q_fossil))
		tree <- read.tree(FBD_temp)
		useBMT=FALSE		
		
		}
		
	if (q_f == 1)  {
		q_fossil <- 20/sum(tree$edge.length)
		temp_tree = sprintf("temp_%s_n.tree", rep_number)
		write.tree(tree, temp_tree)
		FBD_temp = sprintf("FBDtemp_%s_n.tree", rep_number)
		system(paste('python simulate_fossils_on_tree.py ',temp_tree, FBD_temp, q_fossil))
		tree <- read.tree(FBD_temp)
		foss_names<-getExtinct(tree)
		# edges that lead to fossils
		edg_foss <- which(tree$edge[,2] %in% which(tree$tip.label %in% foss_names) )
				
		## based on edges age
	
		edges=cbind(c(1:(n=max(tree$edge[,1])-1)),tree$edge,nodeHeights(tree))
		edges=edges[sort.list(edges[,3]),]
		colnames(edges)<-c("edge", "node1", "node2", "height_n1", "heigh_n2")
		only_foss_edges <- edges[ which(edges[,1] %in% edg_foss) ,] 
		max_fossil <- only_foss_edges [ which(only_foss_edges[,5] == min(only_foss_edges[,5]) ) ,3]
		oldest_fossil <- tree$tip.label [ max_fossil]
		foss_rem <- foss_names [- which(foss_names==oldest_fossil)]
		oldest_fossil_tree <- drop.tip(tree, foss_rem)
		tree <- oldest_fossil_tree
			}
		
			
	if (file.exists(temp_tree)) file.remove(temp_tree)
	if (file.exists(FBD_temp)) file.remove(FBD_temp)
	

	## This line takes the max node height (root) and divide to "scale" trees to a max height of 1. 
	tree$edge.length<- tree$edge.length/max(nodeHeights(tree))
	tree <- ladderize(tree)
	mod_tree <- tree
	ntips = length(tree$tip.label)

	ttree=make.era.map(tree,0)
	phy=tree
	# adjusting it to our data

	for (m in c(1:length(phy$edge.length))){
	    names(ttree$maps[[m]])=m
	    }
	
	test=matrix(0,ncol=length(phy$edge.length), nrow=length(phy$edge.length))
	rownames(test)=rownames(ttree$mapped.edge)
	colnames(test)=c(1:length(phy$edge.length))
	
	diag(test)=ttree$edge.length
	ttree$mapped.edge=test
	
	# verifying is done !
	col = rep("black",length(phy$edge.length))
	names(col)=c(1:length(phy$edge.length))
	#plotSimmap(ttree, colors=col)
	
	# preparing rates for the simualtion 
	sigmas2= rep(s2,length(phy$edge.length))
	mu0s= rep(mu0,length(phy$edge.length))
	
	## SHIFTS IN RATES 
	print("shifts in rates")
	
	rnd_sample_tax = sample((1:ntips)+ntips, ntips)
	chosen_nodes =c()
	r_multipliers = runif(n_shifts,8,16)
	j=1
	if (n_shifts>=1){
		for (s in 1:n_shifts){
			ind_desc_nodes=c()
			while (length(ind_desc_nodes[ind_desc_nodes<ntips])<15 | length(ind_desc_nodes[ind_desc_nodes<ntips])>25 ){
				i = rnd_sample_tax[j]
				ind_desc_nodes = getDescendants(tree, node=i)
				ind_desc_edges = which(tree$edge[,2]%in% ind_desc_nodes == T)
				ind_desc_edges = c(ind_desc_edges, which(tree$edge[,2]==i))
				j=j+1
				#print(c(i, length(intersect(chosen_nodes,ind_desc_nodes))))
				
				if (length(intersect(chosen_nodes,ind_desc_nodes))>0){ind_desc_nodes=c(1)}			
			}
			#print(c(i, length(intersect(chosen_nodes,ind_desc_nodes))))
			chosen_nodes= append(chosen_nodes,ind_desc_nodes)
			sigmas2[ind_desc_edges] = s2*r_multipliers[s] # (xfold*s)
		}
	}else{ind_desc_edges=c(0)}
	
	names(sigmas2)  = NULL # names(col)#  #names_l

	# SHIFTS in TRENDS
	
	print("shifts in trends")
	
	rnd_sample_tax = sample((1:ntips)+ntips, ntips)
	chosen_nodes =c()
	r_multipliers = rnorm(n_shifts_mu0, 2,0.5)*sample(c(-1,1),n_shifts_mu0, replace=T)
	j=1
	if (n_shifts_mu0>=1){
		for (s in 1:n_shifts_mu0){
			ind_desc_nodes=c()
			#while (length(ind_desc_edges)<15 | length(ind_desc_edges)>25 ){
			while (length(ind_desc_nodes[ind_desc_nodes<ntips])<15 | length(ind_desc_nodes[ind_desc_nodes<ntips])>25 ){ 
				i = rnd_sample_tax[j]
				ind_desc_nodes = getDescendants(tree, node=i)
				ind_desc_edges = which(tree$edge[,2]%in% ind_desc_nodes == T)
				ind_desc_edges = c(ind_desc_edges, which(tree$edge[,2]==i))
				j=j+1
				if (length(intersect(chosen_nodes,ind_desc_nodes))>0){ind_desc_nodes=c(1)}			
		
			}
			chosen_nodes= append(chosen_nodes,ind_desc_nodes)
			mu0s[ind_desc_edges] = r_multipliers[s] # (xfold*s)
		}
	}else{ind_desc_edges=c(0)}


	names(mu0s)  = NULL # names(col)#  #names_l
	
	full_data <- simBMT_shifts(ttree, sig2=sigmas2, a=0, nsim=1,mu=mu0s)
	
		
	return(c(tree,full_data,sigmas2,mu0s,original_simulated_tree, q_fossil, ttree))
}


get_time <- function(){
	return(sum(as.numeric(strsplit(format(Sys.time(), "%X"),":")[[1]])*c(3600,60,1)))
}


#################### RUN SIMULATION
#S = sim_data(ntips=ntips,s2=sig2,xfold=xfold,n_shifts=n_shifts)
#S = sim_data_bmt(ntips=ntips,s2=0.1,m0=0.75,root_value=0)

if (use_RDA==0) { 	
	if (treefile != "NA") { 

		print("Reading data")

		t<- read.nexus(treefile)
		if (class(t)=="phylo") { 
			t <- t
			original_simulated_tree <-t } 
			name_tag=""
	
		if (class(t) == "multiPhylo") { 
			t <- t[[tindex]]
			original_simulated_tree <-t
			name_tag=paste("_tree_",tindex, sep="")
			} 

		if (remF == 1) {
			name_tag = paste(name_tag,"_noFossil_",sep="")
			print("Not using BMT and removing fossils")
			useBMT=FALSE 
			fossil_tips <- setdiff(t$tip.label, getExtant(t))
			extant_tree <- drop.tip(t, fossil_tips)
			t<- extant_tree
			} 

		trait <- read.table(dfile, header=F,row.names=1)
		treetrait <- treedata(t,trait) # match trait data and phylogeny
		tree <- treetrait$phy
		tree$edge.length [which(tree$edge.length==0) ]<- 0.00000001
		if (log_trait_data==10){
			data_dataframe <- log10(treetrait$data)
		} else if (log_trait_data==1){
			data_dataframe <- log(treetrait$data)
		} else{
			data_dataframe <- treetrait$data
		}

		data<- as.vector(data_dataframe)
		data <- data * rescale_trait_data
		names(data)<- rownames(data_dataframe)

		ntips=tree$Nnode+1
		D <- build_table(tree, ntips,data)
		BR= branching.times(tree) # sorted from root to most recent
		dist_from_root = max(BR)-BR
		true_sigmas <- rep(1,length(tree$edge.length))
		rda_exist=1 ## attention here is assuming RDA exist but it doesn't when loads the files

	} else {

	worked<-0
	while (worked==0) {
		S <- try(sim_data_bmt_shifts(s2=sig2,mu0=mu0,xfold_sig2=xfold,n_shifts_sig2=n_shifts,n_shifts_mu0=n_shifts_mu0), silent=F)
		if ('try-error' %in% class(S)) worked <-0
		else worked <-1
		} 
	tree= S[[1]]
	ntips=tree$Nnode+1
	full_data= S[[2]]
	data <- S[[2]][names(S[[2]]) %in% tree$tip.label]
	true_anc = S[[2]][names(S[[2]]) %in% tree$edge]
	names_true_anc = names(true_anc)
	D <- build_table(tree, ntips,data)
	BR= branching.times(tree) # sorted from root to most recent
	dist_from_root = max(BR)-BR
	#phenogram(tree, full_data,add=F, col=rep("black",length(full_data)))
	alter_ind= S[[3]]
	true_sigmas = S[[3]]
	true_mu0s = S[[4]]
	original_simulated_tree <- S[[5]]
	q_fossil= S[[6]]
	ttree= S[[7]]
	name_tag=""
	rda_exist=0

	}


} else if (extant_only==1) {
	tree=S[[1]]
	original_simulated_tree <- S[[5]]
	q_fossil= S[[6]]
	ttree= S[[7]]
	fossil_tips <- setdiff(tree$tip.label, getExtant(tree))
	extant_tree <- drop.tip(tree, fossil_tips)
	extant_tips <- getExtant(tree)
	ind_extant <- which(tree$tip.label %in% extant_tips)
	remove_edge <- vector()
	remove_node <- vector()

	for (i in 1:nrow(tree$edge)) { 
		if (length ( intersect(getDescendants(tree, tree$edge[i,2]), ind_extant )) == 0 ) {
			remove_edge <- append(remove_edge, i)
			remove_node <- append(remove_node, tree$edge[i,1])
				} 
	}
	full_data= S[[2]]
	true_sigmas = S[[3]][-remove_edge]
	true_mu0s = S[[4]][-remove_edge]
	true_anc = S[[2]][names(S[[2]]) %in% tree$edge]
	true_anc = true_anc [- which(names(true_anc) %in% unique(remove_node)) ]
	names_true_anc = names(true_anc)
	tree<- extant_tree
	data<- data[-which(names(data) %in% fossil_tips)]
	name_tag="_extant_only"
	ntips=tree$Nnode+1
	D <- build_table(tree, ntips,data)
	BR= branching.times(tree) # sorted from root to most recent
	dist_from_root = max(BR)-BR
}else{
	tree= S[[1]]
	ntips=tree$Nnode+1
	full_data= S[[2]]
	data <- S[[2]][names(S[[2]]) %in% tree$tip.label]
	true_anc = S[[2]][names(S[[2]]) %in% tree$edge]
	names_true_anc = names(true_anc)
	D <- build_table(tree, ntips,data)
	BR= branching.times(tree) # sorted from root to most recent
	dist_from_root = max(BR)-BR
	#phenogram(tree, full_data,add=F, col=rep("black",length(full_data)))
	alter_ind= S[[3]]
	true_sigmas = S[[3]]
	true_mu0s = S[[4]]
	original_simulated_tree <- S[[5]]
	q_fossil= S[[6]]
	ttree= S[[7]]
	name_tag=""
	rda_exist=0
	
}

# define output name
if (treefile != "NA"){
	filename = sprintf("%s%s%s", opt$options$tindex, opt$options$out,name_tag)
} else {
	filename = sprintf("sim_%s%s%s", opt$options$j, opt$options$out,name_tag)
}


### LOG rel err
prior_tbl = get_calibration_tbl(D,root_calibration)


## here create the txt with all the true_ objects

if (rda_exist == 0) { 

l_obj<-length(c(true_sigmas, true_mu0s, true_anc))
res_txt <- matrix(ncol=(l_obj)+4, nrow=1)
TE=as.matrix(tree$edge)
length(TE[,1])
length(true_sigmas)
colnames(res_txt) <- c("n_shifts_sig2", "n_shifts_mu0","xfold","q_fossil", paste("True_sig2", 1:length(TE[,1]),sep="_"),paste("True_mu0", 1:length( TE[,1]),sep="_"), paste("True_anc",D[,1],sep="_"))
res_txt[1,1] <- n_shifts
res_txt[1,2] <- n_shifts_mu0
res_txt[1,3] <- xfold
res_txt[1,4] <- q_fossil
res_txt[1,5:ncol(res_txt)] <- c(true_sigmas, true_mu0s, true_anc)


logfile4 = sprintf("%s_truestates.txt", filename)
write.table(res_txt, file=paste(logfile4), row.names=FALSE, quote=FALSE, sep="\t")
}



## save the whole R image 
res_all_files = sprintf("%s.rda", filename)
save.image(res_all_files)


print("Starting MCMC...")
t1 = get_time()

logfile3 = sprintf("%s.log", filename)

mcmc.gibbs4(tree, data, D, prior_tbl, true_rate=true_sigmas, ngen= nGibbs_gen,bdmcmc_freq=0.75,logfile=logfile3,update_sig_freq=0.5,
	sample=Gibbs_sample,print_freq=print_f,dynamicPlot=dynamicPlot,useTrend=useBMT) 
cat("\nTime Gibbs:", get_time() - t1, "\n", sep="\t")


# plot res if specified 
resfile = sprintf("%s.pdf", filename)

if (plot_res==1 && treefile == "NA") { 
	plot_results <- function(filename,plot_phenogram=T,rdafile=0,logfile=0){
		require(phytools)
		library(geiger)
		library(TreeSim)
		library(scales)
		# load original data
		if (rdafile==0){
			robj_file = sprintf("%s.rda", filename)
		}else{
			robj_file = rdafile
		}
	
		load(robj_file)

		# load MCMC samples
		if (logfile==0){
			logfile3= sprintf("%s.log", filename)
		}else{
			logfile3=logfile
		}
	
	
		pdf(file=resfile,width=25*0.75, height=15*0.75)
		par(mfrow=c(2,4))

		col = rep("black",length(tree$edge.length))
		names(col)=c(1:length(tree$edge.length))
		plot.phylo(tree, edge.width=true_sigmas*2,show.tip.label = F, main="True rates")
		edge_w=abs(full_data[ttree$edge[,1]] - full_data[ttree$edge[,2]])/sqrt(tree$edge.length)
		#plot.phylo(tree, edge.width=edge_w*3,show.tip.label = F,main="Empirical rates")
		#phenogram(ttree, full_data,add=F, col=col, main="Traigram")

		out_tbl = read.table(logfile3,header=T)
		ind_sig2_col = grep('sig2_', colnames(out_tbl), value=F)
		burnin= round(0.5*dim(out_tbl)[1])
		rates_temp_median = apply(out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=median,2)
		rates_temp_mean = apply(out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=mean,2)
		MAPE_median = mean(abs((rates_temp_median-true_sigmas))/true_sigmas)
		sdAPE_median = sd(abs((rates_temp_median-true_sigmas)))
		MAPE_mean = mean(abs((rates_temp_mean-true_sigmas))/true_sigmas)
		sdAPE_mean = sd(abs((rates_temp_mean-true_sigmas)))
		estK_all =out_tbl$K_sig2[burnin:dim(out_tbl)[1]]
		h=hist(estK_all, breaks=0:max(estK_all)+0.5,plot=F)
		best_K = which(h$counts==max(h$counts))
		estK = round(quantile(estK_all,probs = c(0.025,0.975)))
		rates_min = apply(X=out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=quantile, probs=0.025,2)
		rates_max = apply(X=out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=quantile, probs=0.975,2)
		plot.phylo(tree, edge.width=rates_temp_mean*2, main=paste("Estimated rates,",sprintf("K: %s (%s-%s)",best_K,estK[1],estK[2])),show.tip.label = F)
		C = h$counts
		names(C) = 1:max(estK_all)
		barplot(C,main="Estimated number of shifts - Sig2")
		plot(true_sigmas,  main=paste("MAPE_avg:",round(MAPE_mean,3),"MAPE_med:",round(MAPE_median,3)),pch=19,col="darkblue",
		ylim=c(0,max(c(true_sigmas,rates_max))),type="n",lwd=3)
		segments(x0= 1:length(rates_min), y0=rates_min, x1=1:length(rates_min), y1=rates_max,col=alpha("darkred",0.4),lwd=2)
		points(rates_temp_mean,pch=20,col="darkred")
		points(true_sigmas, col="darkblue",type="l",lwd=3)
	

		ind_anc_col = grep('anc_', colnames(out_tbl), value=F)
		anc_mean = apply(out_tbl[burnin:dim(out_tbl)[1], ind_anc_col],FUN=mean,2)
		mape = round(mean(abs(anc_mean-true_anc)/(max(true_anc)-min(true_anc))),3)

		estKmu0_all =out_tbl$K_mu0[burnin:dim(out_tbl)[1]]
		h=hist(estKmu0_all, breaks=0:max(estKmu0_all)+0.5,plot=F)
		best_K = which(h$counts==max(h$counts))
		estKmu0 = round(quantile(estKmu0_all,probs = c(0.025,0.975)))
		ind_mu0_col = grep('mu0_', colnames(out_tbl), value=F)
		mu0_temp_mean = (apply(out_tbl[burnin:dim(out_tbl)[1], ind_mu0_col],FUN=mean,2))
		mu0_min = apply(X=out_tbl[burnin:dim(out_tbl)[1], ind_mu0_col],FUN=quantile, probs=0.025,2)
		mu0_max = apply(X=out_tbl[burnin:dim(out_tbl)[1], ind_mu0_col],FUN=quantile, probs=0.975,2)
		#hist(mu0_temp_median, main=paste("Estimated mu0s",sprintf("(K: %s-%s)",estKmu0[1],estKmu0[2])))
		plot(true_mu0s,  main=paste("mu0s, ", sprintf("K: %s (%s-%s)",best_K,estKmu0[1],estKmu0[2])),pch=21,col="darkblue",bg="blue",
		ylim=c(min(c(true_mu0s,mu0_min))-abs(0.5*min(c(true_mu0s,mu0_min))),max(c(true_mu0s,mu0_max))+abs(0.5*max(c(true_mu0s,mu0_max)))))
		points(mu0_temp_mean,pch=21,col="darkred",bg="red")
		segments(x0= 1:length(mu0_min), y0=mu0_min, x1=1:length(mu0_min), y1=mu0_max,col=alpha("darkred",0.4),lwd=2)
		points(mu0_temp_mean,pch=20,col="darkred")
		points(true_mu0s, col="darkblue",type="l",lwd=3)
		names(anc_mean) = names(true_anc)
		est_anc_states = c(data,anc_mean)
		C = h$counts
		names(C) = 1:max(estKmu0_all)
		barplot(C,main="Estimated number of shifts - Mu0")
		if (plot_phenogram==T){
			phenogram(tree, full_data,add=F, col=alpha("darkblue",0.2), main="Traigram",ylim=c(min(est_anc_states,full_data),max(est_anc_states,full_data)), spread.labels=F)
			phenogram(tree, est_anc_states,add=T, col="darkred", main="Traigram", spread.labels=F)		
		}
		R2_anc = summary(lm(anc_mean ~ true_anc))$r.squared
		plot(anc_mean, true_anc, xlab="estimated anc states", ylab="true anc states",main=paste("Anc states - R2: ", round(R2_anc,3)))	
		abline(coef=c(0,1),lty=2)
	
		R2_sig2 = mean(abs(rates_temp_mean - true_sigmas))
		R2_mu0 = mean(abs(mu0_temp_mean - true_mu0s))
		print(R2_sig2)
		print(R2_mu0)
		n<-dev.off()
	
	
	}
	plot_results(filename)
}


if (plot_res==1 && treefile != "NA"){
	## Modified traitgram 
	phenogram_invTime<-function(tree,x,fsize=1.0,ftype="reg",colors=NULL,axes=list(),add=FALSE,...){
		## get optional arguments
		if(hasArg(xlim)) xlim<-list(...)$xlim
		else xlim<-NULL
		if(hasArg(ylim)) ylim<-list(...)$ylim
		else ylim<-NULL
		if(hasArg(log)) log<-list(...)$log
		else log<-""
		if(hasArg(main)) main<-list(...)$main
		else main<-NULL
		if(hasArg(sub)) sub<-list(...)$sub
		else sub<-NULL
		if(hasArg(xlab)) xlab<-list(...)$xlab
		else xlab<-"time"
		if(hasArg(ylab)) ylab<-list(...)$ylab
		else ylab<-"phenotype"
		if(hasArg(asp)) asp<-list(...)$asp
		else asp<-NA
		if(hasArg(type)) type<-list(...)$type
		else type<-"l"
		if(hasArg(lty)) lty<-list(...)$lty
		else lty<-1
		if(hasArg(lwd)) lwd<-list(...)$lwd
		else lwd<-2
		if(hasArg(offset)) offset<-list(...)$offset
		else offset<-0.2
		if(hasArg(offsetFudge)) offsetFudge<-list(...)$offsetFudge
		else offsetFudge<-1.37
		if(hasArg(digits)) digits<-list(...)$digits
		else digits<-2
		if(hasArg(nticks)) nticks<-list(...)$nticks
		else nticks<-5
		if(hasArg(spread.labels)) spread.labels<-list(...)$spread.labels
		else spread.labels<-TRUE
		if(ftype=="off") spread.labels<-FALSE
		if(hasArg(spread.cost)) spread.cost<-list(...)$spread.cost
		else spread.cost<-c(1,0.4)
		if(hasArg(spread.range)) spread.range<-list(...)$spread.range
		else spread.range<-range(x)
		if(hasArg(link)) link<-list(...)$link
		else link<-if(spread.labels) 0.1*max(nodeHeights(tree)) else 0
		if(hasArg(hold)) hold<-list(...)$hold
		else hold<-TRUE
		if(hasArg(quiet)) quiet<-list(...)$quiet
		else quiet<-FALSE
		## end optional arguments
		# check tree
		if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
		# check font
		ftype<-which(c("off","reg","b","i","bi")==ftype)-1
		if(!ftype&&!add) fsize=0 
		H<-nodeHeights(tree)
		if(length(x)<(length(tree$tip)+tree$Nnode))
			x<-c(x,fastAnc(tree,x))
		else
			x<-c(x[tree$tip.label],x[as.character(length(tree$tip)+1:tree$Nnode)])
		x[1:length(tree$tip)]<-x[tree$tip.label]
		names(x)[1:length(tree$tip)]<-1:length(tree$tip)
		X<-matrix(x[as.character(tree$edge)],nrow(tree$edge),ncol(tree$edge))
		# legacy 'axes' argument trumps ylim & xlim from optional (...)
		if(is.null(axes$trait)&&is.null(ylim)) ylim<-c(min(x),max(x))
		else if(!is.null(axes$trait)) ylim<-axes$trait
		if(!is.null(axes$time)) xlim<-axes$time
		if(!add&&is.null(xlim)){
			pp<-par("pin")[1]
			sw<-fsize*(max(strwidth(tree$tip.label,units="inches")))+offsetFudge*offset*fsize*strwidth("W",units="inches")
			alp<-optimize(function(a,H,link,sw,pp) (a*1.04*(max(H)+link)+sw-pp)^2,H=H,link=link,sw=sw,pp=pp,interval=c(0,1e6))$minimum
			xlim<-c(min(H),max(H)+link+sw/alp)
		}
		if(!quiet&&Ntip(tree)>=40&&spread.labels){ 
			cat("Optimizing the positions of the tip labels...\n")
			flush.console()
		}
		## matrix for tip coordinates
		tip.coords<-matrix(NA,Ntip(tree),2,dimnames=list(tree$tip.label,c("x","y")))
		if(hold) null<-dev.hold()
		if(is.null(tree$maps)){
			if(is.null(colors)) colors<-"black"
			if(!add){ 
				plot(H[1,],X[1,],type=type,lwd=lwd,lty=lty,col=colors,xlim=xlim,ylim=ylim,log=log,asp=asp,xlab="",ylab="",frame=FALSE, axes=FALSE)
				if(spread.labels) tt<-spreadlabels(tree,x,fsize=fsize,cost=spread.cost,range=spread.range) else tt<-x[1:length(tree$tip)]
				if(tree$edge[1,2]<=length(tree$tip)){
					if(fsize&&!add){
						text(tree$tip.label[tree$edge[1,2]],x=H[1,2]+link,y=tt[tree$edge[1,2]],cex=fsize,font=ftype,pos=4,offset=offset)
						tip.coords[tree$tip.label[tree$edge[1,2]],]<-c(H[1,2]+link,tt[tree$edge[1,2]])
						if(link>0) lines(x=c(H[1,2],H[1,2]+link),y=c(X[1,2],tt[tree$edge[1,2]]),lty=3)
					}
				}
				s<-2
			} else s<-1
			for(i in s:nrow(H)){ 
				lines(H[i,],X[i,],type=type,lwd=lwd,lty=lty,col=colors)
				if(tree$edge[i,2]<=length(tree$tip)){
					if(fsize&&!add){ 
						text(tree$tip.label[tree$edge[i,2]],x=H[i,2]+link,y=tt[tree$edge[i,2]],cex=fsize,font=ftype,pos=4,offset=offset)
						tip.coords[tree$tip.label[tree$edge[i,2]],]<-c(H[i,2]+link,tt[tree$edge[i,2]])
						if(link>0) lines(x=c(H[i,2],H[i,2]+link),y=c(X[i,2],tt[tree$edge[i,2]]),lty=3)
					}
				}
			}
		} else {
			if(is.null(colors)){
				nn<-sort(unique(c(getStates(tree,"tips"),getStates(tree,"nodes"))))
				colors<-setNames(palette()[1:length(nn)],nn)
			}
			for(i in 1:nrow(H)){
				y<-H[i,1]
				m<-diff(X[i,])/diff(H[i,])
				for(j in 1:length(tree$maps[[i]])){
					a<-c(y,y+tree$maps[[i]][j])
					b<-m*(a-H[i,1])+X[i,1]
					if(i==1&&j==1&&!add) {
						plot(a,b,col=colors[names(tree$maps[[i]])[j]],type=type,lwd=lwd,lty=lty,xlim=xlim,ylim=ylim,log=log,asp=asp,axes=FALSE,xlab="",ylab="")
						if(spread.labels) tt<-spreadlabels(tree,x[1:length(tree$tip)],fsize=fsize,cost=spread.cost,range=spread.range) else tt<-x[1:length(tree$tip)]
					} else lines(a,b,col=colors[names(tree$maps[[i]])[j]],lwd=lwd,lty=lty,type=type)
					y<-a[2]
				}
				if(tree$edge[i,2]<=length(tree$tip)){
					if(fsize&&!add){ 
						text(tree$tip.label[tree$edge[i,2]],x=H[i,2]+link,y=tt[tree$edge[i,2]],cex=fsize,font=ftype,pos=4,offset=offset)
						tip.coords[tree$tip.label[tree$edge[i,2]],]<-c(H[i,2]+link,tt[tree$edge[i,2]])
						if(link>0) lines(x=c(H[i,2],H[i,2]+link),y=c(X[i,2],tt[tree$edge[i,2]]),lty=3)
					}
				}
			}
		}
		if(!add){
			at<-round(0:(nticks-1)*max(H)/(nticks-1),digits)
			lab_at<- rev(at)
			axis(1,at=at, labels=lab_at); axis(2); title(xlab=xlab,ylab=ylab,main=main,sub=sub)
		}
		if(hold) null<-dev.flush()
		xx<-setNames(c(H[1,1],H[,2]),c(tree$edge[1,1],tree$edge[,2]))
		xx<-xx[order(as.numeric(names(xx)))]
		yy<-setNames(c(X[1,1],X[,2]),c(tree$edge[1,1],tree$edge[,2]))
		yy<-yy[order(as.numeric(names(yy)))]
		PP<-list(type="phenogram",use.edge.length=TRUE,node.pos=1,
			show.tip.label=if(ftype!="off") TRUE else FALSE,show.node.label=FALSE,
			font=ftype,cex=fsize,adj=0,srt=NULL,no.margin=FALSE,label.offset=offset,
			x.lim=par()$usr[1:2],y.lim=par()$usr[3:4],
			direction=NULL,tip.color="black",Ntip=Ntip(tree),Nnode=tree$Nnode,
			edge=tree$edge,xx=xx,yy=yy)
		assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
		invisible(tip.coords)
	}

	## function to spread labels
	## written by Liam J. Revell 2013, 2014, 2016
	spreadlabels<-function(tree,x,fsize=1,cost=c(1,1),range=NULL){
		if(is.null(range)) range<-range(x)
		yy<-x[1:Ntip(tree)]
		zz<-setNames((rank(yy,ties.method="random")-1)/(length(yy)-1)*diff(range(yy))+range(yy)[1],names(yy))
		mm<-max(fsize*strheight(tree$tip.label))
		ff<-function(zz,yy,cost,mo=1,ms=1){
			ZZ<-cbind(zz-mm/2,zz+mm/2)
			ZZ<-ZZ[order(zz),]
			oo<-0
			for(i in 2:nrow(ZZ)) 
				oo<-if(ZZ[i-1,2]>ZZ[i,1]) oo<-oo+ZZ[i-1,2]-ZZ[i,1] else oo<-oo
			pp<-sum((zz-yy)^2)
			oo<-if(oo<(1e-6*diff(par()$usr[3:4]))) 0 else oo
			pp<-if(pp<(1e-6*diff(par()$usr[3:4]))) 0 else pp
			oo/mo*cost[1]+pp/ms*cost[2]
		}
		mo<-ff(yy,zz,cost=c(1,0))
		ms<-ff(yy,zz,cost=c(0,1))
		if(mo==0&&ms==0) return(yy)
		else {
			rr<-optim(zz,ff,yy=yy,mo=mo,ms=ms,cost=cost,method="L-BFGS-B",lower=rep(range[1],length(yy)),upper=rep(range[2],length(yy)))
	return(rr$par)
		}
	}



	plot_results <- function(fRDA, fLOG, resfile , plot_phenogram=T){	
		require(phytools)
		#require(diversitree)
		require(methods)
		library(geiger)
		library(TreeSim)
		library(scales)
		library(plotrix)
		load(fRDA)
			
		out_tbl = read.table(fLOG,header=T)
		ind_sig2_col = grep('sig2_', colnames(out_tbl), value=F)
		burnin= round(0.25*dim(out_tbl)[1])

		pdf(file=resfile,width=15*0.75, height=25*0.75)
	
		# Rates section
	
		rates_temp_median = apply(out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=median,2)
		rates_temp_mean = apply(out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=mean,2)
		estK_all =out_tbl$K_sig2[burnin:dim(out_tbl)[1]]
		h=hist(estK_all, breaks=0:max(estK_all)+0.5,plot=F)
		best_K = which(h$counts==max(h$counts))
		estK = round(quantile(estK_all,probs = c(0.025,0.975)))
		rates_min = apply(X=out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=quantile, probs=0.025,2)
		rates_max = apply(X=out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=quantile, probs=0.975,2)

		par(mfcol=c(3,2))

		temp<-log(rates_temp_mean) + max(abs(log(rates_temp_mean))) + 0.5

		# Rates are log-transformed to avoid extremes, then make sure that values are positive, then 
		# add an arbitrary small number to avoid zeros for the cex - edge width parameter in plotting.
	
		col = rep("black",length(tree$edge.length))
		names(col)=c(1:length(tree$edge.length))
		plot.phylo(tree, edge.width=temp, main=paste("Estimated rates,",sprintf("K: %s (%s-%s)",best_K,estK[1],estK[2])),show.tip.label = T, cex=0.6, align.tip.label=T)
	
		C = h$counts
		names(C) = 1:max(estK_all)
		barplot(C,main="Estimated number of shifts - Sig2")
	
		plot(rates_temp_mean,  main=paste("sig2s, ", sprintf("K: %s (%s-%s)",best_K,estK[1],estK[2])),pch=21,col="darkblue",bg="blue", ylim=c(min(rates_min)*0.5,max(rates_max)*1.5),type="n")
		points(rates_temp_mean,pch=21,col="darkred",bg="red")
		segments(x0= 1:length(rates_min), y0=rates_min, x1=1:length(rates_min), y1=rates_max,col=alpha("darkred",0.4),lwd=2)
	
	
		# Trends section
	
		estKmu0_all =out_tbl$K_mu0[burnin:dim(out_tbl)[1]]
		h=hist(estKmu0_all, breaks=0:max(estKmu0_all)+0.5,plot=F)
		best_K = which(h$counts==max(h$counts))
		estKmu0 = round(quantile(estKmu0_all,probs = c(0.025,0.975)))
	
		ind_mu0_col = grep('mu0_', colnames(out_tbl), value=F)
		mu0_temp_mean = (apply(out_tbl[burnin:dim(out_tbl)[1], ind_mu0_col],FUN=mean,2))
		mu0_min = apply(X=out_tbl[burnin:dim(out_tbl)[1], ind_mu0_col],FUN=quantile, probs=0.025,2)
		mu0_max = apply(X=out_tbl[burnin:dim(out_tbl)[1], ind_mu0_col],FUN=quantile, probs=0.975,2)
	
		rbPal <- colorRampPalette(c("#2166ac","gray","#b2182b"))
		# maximum color scale for trend 
		xa <- max(abs(c (mu0_temp_mean)))+ (0.1* max(abs(c (mu0_temp_mean))))
		min_col <- -xa
		max_col <- xa
		beta_shape <- 0.5
		test<- qbeta(seq(0,1, by=0.05), beta_shape, beta_shape, lower.tail = TRUE, log.p = FALSE)
		x = test * (max_col-min_col) 
		breaksList = x + min_col
		
		edge_cols <- rbPal(length(breaksList))[as.numeric(cut(sort(mu0_temp_mean),breaks = breaksList))]
		names(edge_cols) <- order(mu0_temp_mean)
		edge_cols <- edge_cols [order(as.numeric(names(edge_cols)))]
		names(edge_cols)<- NULL
		plot.phylo(tree, edge.width=2, main=paste("Estimated trends,",sprintf("K: %s (%s-%s)",best_K,estK[1],estK[2])),show.tip.label = F, edge.color=edge_cols)
		testcol<-rbPal(n=length(breaksList))
	
		col.labels<-c(format(min_col, digits = 4), "","0","", format(max_col, digits = 4))
	  	limit_plot<-max(axisPhylo())
	  	color.legend(0,0.15,0+(limit_plot)/10,25,col.labels,rect.col=testcol, gradient="y", cex=0.8, align="lt")

		C = h$counts
		names(C) = 1:max(estKmu0_all)
		barplot(C,main="Estimated number of shifts - Mu0")
		
		plot(mu0_temp_mean,  main=paste("mu0s, ", sprintf("K: %s (%s-%s)",best_K,estKmu0[1],estKmu0[2])),pch=21,col="darkblue",bg="blue", ylim=c(min(mu0_min)-abs(0.5*min(mu0_min)),max(mu0_max)+abs(0.5*min(mu0_min))),type="n")
		points(mu0_temp_mean,pch=21,col="darkred",bg="red")
		segments(x0= 1:length(mu0_min), y0=mu0_min, x1=1:length(mu0_min), y1=mu0_max,col=alpha("darkred",0.4),lwd=2)
		
		## Ancestral states section
		par(mfcol=c(2,1))

		ind_anc_col = grep('anc_', colnames(out_tbl), value=F)
		anc_mean = apply(out_tbl[burnin:dim(out_tbl)[1], ind_anc_col],FUN=mean,2)

		names(anc_mean) = seq(from=ntips+1, to=max(tree$edge[,2]), by=1)
		est_anc_states = c(data,anc_mean)
	
		use_mcmc <- out_tbl[burnin:dim(out_tbl)[1], ind_anc_col]
		sample_size <- min(100, dim(use_mcmc)[1])
		sample_mcmc <- use_mcmc [ sample(1:dim(use_mcmc)[1], size=sample_size, replace = FALSE) , ]
		colnames(sample_mcmc) = seq(from=ntips+1, to=max(tree$edge[,2]), by=1)
		sample_mcmc <- as.matrix(sample_mcmc)
	
		if (plot_phenogram==T){
	
		# Comming back to the original values of the trait, and applying the same transform/scaling used for MCMC
	
		if(log_trait_data == 0) {
		   sample_mcmc_raw <- sample_mcmc
		   est_anc_states_mean <- est_anc_states
		   data_raw <- data
		  	}
		if(log_trait_data == 1){ 
		  sample_mcmc_raw <- exp(sample_mcmc)	 
		  est_anc_states_mean <- exp(est_anc_states)
		  data_raw <- exp(data)	
		   }
		if(log_trait_data == 10){
		  sample_mcmc_raw <- 10^sample_mcmc
		  est_anc_states_mean <- 10^est_anc_states
		  data_raw <- 10^data	
	    		}		
		data_raw = data/rescale_trait_data
		sample_mcmc_raw = sample_mcmc/rescale_trait_data
		est_anc_states_mean <- est_anc_states/rescale_trait_data
	 
		# mean values phenogram 
	
		color=adjustcolor( "darkred", alpha.f = 1)
		phenogram_invTime(tree, est_anc_states_mean, col=color, main="Traigram", spread.labels=F)
	 
		# mcmc sample of phenograms
	
		alphacol = 1/sample_size + sample_size*0.001
	
		color=adjustcolor( "darkred", alpha.f = alphacol)
		temp = c(data_raw, sample_mcmc_raw[1,])	
		phenogram_invTime(tree, temp, col=color, main="Traigram", spread.labels=F, ftype="off", ylim=c(min(sample_mcmc_raw), max(sample_mcmc_raw)))
	
			for (s in 2:dim(sample_mcmc)[1]) { 
			temp = c(data_raw, sample_mcmc_raw[s,])	
			phenogram_invTime(tree, temp, col=color, main="Traigram", spread.labels=F, ftype="off", add=T)
			}

		n<-dev.off()
	
		}

	}
	
	plot_results(res_all_files, logfile3, resfile)

	
}