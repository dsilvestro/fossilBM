library(scales)
require(phytools)
library(geiger)
library(adephylo)



read_and_transform_data <- function(treefile="", datafile="", tree_obj="", data_obj="",
									rm_extinct=FALSE, tindex=1, drop_na=FALSE,
									log_trait_data=0, rescale_trait_data=1,root_calibration=c(0,100),
									partition_file="",zero_br=0, dist_from_the_root=c(),
                                    state_tbl=c(), z_tranform=FALSE){
		
		if (treefile != ""){
			t<- read.nexus(treefile)
		}
		
		if (class(tree_obj)=="phylo") { 
			t <- tree_obj
			t$edge.length <- t$edge.length + zero_br
			original_simulated_tree <-t 
			name_tag=""
			} 

		if (class(tree_obj) == "multiPhylo") { 
			t <- tree_obj[[tindex]]
			original_simulated_tree <-t
			t$edge.length <- t$edge.length + zero_br
			name_tag=paste("_tree_",tindex, sep="")
			} 
        #
        # if (tree_obj[1] != ""){
        #     t = tree_obj
        #     t$edge.length <- t$edge.length + zero_br
        #     name_tag=""
        # }

		if (rm_extinct) {
			name_tag = paste(name_tag,"_noFossil_",sep="")
			print("Not using BMT and removing fossils")
			useBMT=FALSE 
			fossil_tips <- setdiff(t$tip.label, getExtant(t))
			extant_tree <- drop.tip(t, fossil_tips)
			t<- extant_tree
			} 

		fbm_obj = NULL
		
		if (datafile != ""){
			trait <- read.table(datafile, header=F,row.names=1)
		}else{
			trait = data_obj
		}
		
		if (drop_na){
			trait <- na.omit(trait)
		}
		treetrait <- treedata(t,trait) # match trait data and phylogeny
		tree <- treetrait$phy
		tree$edge.length [which(tree$edge.length==0) ]<- 0.00000001
		fbm_obj$tree <- tree
		if (log_trait_data==10){
			data_dataframe <- log10(treetrait$data)
		} else if (log_trait_data==1){
			data_dataframe <- log(treetrait$data)
		} else{
			data_dataframe <- treetrait$data
		}
		
		data <- as.vector(data_dataframe)
		data <- data * rescale_trait_data
		names(data)<- rownames(data_dataframe)
		fbm_obj$data <- data

		fbm_obj$ntips <- tree$Nnode+1
		fbm_obj$D <- build_table(tree, fbm_obj$ntips,fbm_obj$data)
		BR <- branching.times(tree) # sorted from root to most recent
        if (length(dist_from_the_root) == 0){
            fbm_obj$dist_from_the_root <- c(distRoot(fbm_obj$tree,tips="all",method =c( "patristic")), max(BR)-BR)
        }else{
            fbm_obj$dist_from_the_root = dist_from_the_root
        }
		fbm_obj$dist_from_midpoint <- fbm_obj$dist_from_the_root - 0.5*max(fbm_obj$dist_from_the_root)
		fbm_obj$prior_tbl <- get_calibration_tbl(fbm_obj$D,root_calibration)
		fbm_obj$PartitionFile <- partition_file
		fbm_obj$trait_rescaling <- rescale_trait_data
        fbm_obj$StateTbl = state_tbl
		
        if (z_tranform){
            fbm_obj$original_data <- fbm_obj$data
            fbm_obj$z_sd <- sd(fbm_obj$data, na.rm=T)
            fbm_obj$data <- fbm_obj$data / fbm_obj$z_sd 
            fbm_obj$z_mean <- mean(fbm_obj$data, na.rm=T)
            fbm_obj$data <- fbm_obj$data - fbm_obj$z_mean 
        }
        
		return(fbm_obj)	
}

drop_extinct <- function(tree, rho=0.25){
	dropRandomTaxa = 0.25
	ex_tips <- getExtinct(tree)
	rm_taxa= sample(ex_tips,size=round((1-rho) * length(ex_tips)))
	tree <- drop.tip(tree, rm_taxa)
	return(tree)
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

newlnLike <- function(fbm_obj, vector_tip_root_nodes_values, sigma2,mu0, a0) {
	tree <- fbm_obj$tree
	ntips <- fbm_obj$ntips
	D <- fbm_obj$D
	root_dist <- as.numeric(fbm_obj$dist_from_midpoint)
	names(root_dist) <- names(vector_tip_root_nodes_values)
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
	L1 = dnorm(a_v, mean=(anc_v+mu0[a_ind_v]*vpa_v + a0[a_ind_v]*vpa_v*root_dist[a_ind_v]), sd= sqrt(vpa_v*sigma2[a_ind_v]),log=T)
	L2 = dnorm(b_v, mean=(anc_v+mu0[b_ind_v]*vpb_v + a0[b_ind_v]*vpb_v*root_dist[b_ind_v]), sd= sqrt(vpb_v*sigma2[b_ind_v]),log=T)  
	L = L1+L2
	return(L) 
}

phylo_imputation <-function(tree, vector_tip_root_nodes_values, sigma2,D,mu0,a0,ind_NAtaxa_in_D,anc_of_NAtaxa,brl_NAtaxa_in_D,root_dist,get_expected=0) {
	vec_values = vector_tip_root_nodes_values
	anc_v <- vec_values[anc_of_NAtaxa] # anc value
	if (get_expected==0){
		new_val = rnorm(length(brl_NAtaxa_in_D), mean=(anc_v+mu0[ind_NAtaxa_in_D]*brl_NAtaxa_in_D + a0[ind_NAtaxa_in_D]*root_dist[ind_NAtaxa_in_D]*brl_NAtaxa_in_D), sd= sqrt(brl_NAtaxa_in_D*sigma2[ind_NAtaxa_in_D]))
	}else{
		new_val = anc_v + mu0[ind_NAtaxa_in_D]*brl_NAtaxa_in_D + a0[ind_NAtaxa_in_D]*root_dist[ind_NAtaxa_in_D]*brl_NAtaxa_in_D
	}
	# print(c("mean", (anc_v+mu0[ind_NAtaxa_in_D]*brl_NAtaxa_in_D), "std",sqrt(brl_NAtaxa_in_D*sigma2[ind_NAtaxa_in_D]),new_val))
	# print(new_val)
	return(new_val) 
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

get_mu_t <- function(t, x, mu0, b, a0=0){
	return( a0*t*x + mu0*x +b )
}


runGibbs <- function(fbm_obj,sigma2, vector_tip_root_nodes_values,mu0, a0,get_expected=0) {
	ntips      <- fbm_obj$ntips
	D          <- fbm_obj$D
	prior_tbl  <- fbm_obj$prior_tbl
	#prior_tbl = calc_prior_dist_root_calibration(D,sigma2)
	root_dist <- fbm_obj$dist_from_midpoint
	
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

		#desc_lik = get_joint_mu_s2((a-vpa*mu0[sig_a_ind]), (vpa)*sigma2[sig_a_ind], (b-vpb*mu0[sig_b_ind]), (vpb)*sigma2[sig_b_ind])
		
		m1 <- a - (vpa*mu0[sig_a_ind] + a0[sig_a_ind]*vpa*root_dist[a_ind])
		s1 <- (vpa)*sigma2[sig_a_ind]
		m2 <- b - (vpb*mu0[sig_b_ind] + a0[sig_b_ind]*vpb*root_dist[b_ind])
		s2 <- (vpb)*sigma2[sig_b_ind]
		
		desc_lik = get_joint_mu_s2(m1, s1, m2,  s2)
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
			
			stem_lik_m1 <- stem_val + stem_brl*mu0[sig_stem_ind] + a0[sig_stem_ind]*stem_brl*root_dist[anc_ind]
			stem_lik_s1 <- stem_brl * sigma2[sig_stem_ind]
			
			lik_val = get_joint_mu_s2(stem_lik_m1,stem_lik_s1,desc_lik[1],desc_lik[2])

			post_prm= lik_val #get_joint_mu_s2(lik_val[1],lik_val[2],prior_prm[1], prior_prm[2])
			if (get_expected==1){
				vec_values[anc_ind] = post_prm[1]
			}else{
				vec_values[anc_ind] = rnorm(1, mean=post_prm[1], sd=sqrt(post_prm[2]))
			}
		}
	}
	return(vec_values[D[,1]] )
}

calc_prior <- function(sig2, a, y, mu0, delta_a0, prior_tbl, rate_sig2, sd_mu0) {
	prior_sig2 = sum(dexp(sig2, rate_sig2, log = TRUE) ) #sum(dgamma(sig2, 2,1, log = TRUE) )
	prior_root = sum(dnorm(c(a), mean = prior_tbl[,1], sd = prior_tbl[,2], log = T))
	prior_mu0  = sum(dnorm(mu0[mu0 != 0], mean = 0, sd = sd_mu0, log = T)) # ignore mu0 constrained to 0
	prior_da0  = sum(dnorm(delta_a0, mean = 0, sd = 1, log = T))
	return(prior_sig2+prior_root+prior_mu0+prior_da0)
}

set_model_partitions <- function(fbm_obj,ind_sig2,ind_mu0){
    # ind_sig2 <- rep(1, fbm_obj$tree$Nnode * 2 - 1)
    # ind_mu0  <- rep(1, fbm_obj$tree$Nnode * 2 - 1)
	tbl = read.table(fbm_obj$PartitionFile,h=F,stringsAsFactors=F,fill=T)
	for (i in 1:dim(tbl)[1]){
		tx_tmp = as.vector(unlist(tbl[i,]))
		tx = tx_tmp[tx_tmp != ""]
		mrca = getMRCA(fbm_obj$tree, tx)
		desc = getDescendants(fbm_obj$tree, mrca)
		desc[desc > fbm_obj$ntips] = desc[desc > fbm_obj$ntips] - 1
		ind_sig2[desc] = i+1
		ind_mu0[desc] = i+1
	}
	
	sig2 <- rep(0.05, 1+dim(tbl)[1])
	mu0  <- rep(0,   1+dim(tbl)[1])
	a0  <-  rep(0,   1+dim(tbl)[1])
	return( list(sig2, mu0, a0, ind_sig2, ind_mu0) )
}


set_model_trait_partitions <- function(fbm_obj){
	ind_mu0  <- fbm_obj$StateTbl[,2]
    
    
	mu0  <- rep(0, length(unique(ind_mu0)))
	a0  <-  rep(0, length(unique(ind_mu0)))
	return( list(mu0, a0, ind_mu0) )
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

calc_deathRates <- function (fbm_obj, trait_states, sig2, ind_sig, mu0, ind_mu0, a0,IND_edge,desc_index_list,death_sig2=1){
	tree <- fbm_obj$tree
	D <- fbm_obj$D
	ntips <- fbm_obj$ntips
	list_deathrates=c(0) # can't kill first
	y_temp = runGibbs(fbm_obj,sig2[ind_sig], trait_states,mu0[ind_mu0],a0[ind_mu0],get_expected=1)
	y_temp_FIRST = y_temp
	a= y_temp[1]
	y= y_temp[-1]
	trait_states_ = c(trait_states[1:length(tree$tip.label)], a, y)
	current_lik = sum(newlnLike(fbm_obj, trait_states_, sig2[ind_sig],mu0[ind_mu0],a0[ind_mu0]))
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
			a0_temp = a0[ind_mu0]
		}else{
			new_ind_mu0 = ind_mu0
			new_ind_mu0[new_ind_mu0==i] = ind_of_ancestor		
			rates_temp= sig2[ind_sig]
			mu0_temp = mu0[new_ind_mu0]
			a0_temp = a0[new_ind_mu0]
		}
		y_temp = runGibbs(fbm_obj,rates_temp, trait_states,mu0_temp,a0_temp,get_expected=1)
		a= y_temp[1]
		y= y_temp[-1]
		trait_states_temp = c(trait_states[1:length(tree$tip.label)], a, y)
		list_deathrates[i] = sum(newlnLike(fbm_obj, trait_states_temp, rates_temp,mu0_temp,a0_temp))
		
	}
	list_deathrates = exp(list_deathrates - current_lik)	
	#print(c(current_lik, list_deathrates))
	list_deathrates[1] = 0
	return(list_deathrates)
}

run_BDMCMC <- function(fbm_obj, trait_states, sig2, ind_sig, mu0, ind_mu0, a0 ,IND_edge,new_ind,desc_index_list,death_sig2=1){
	tree <- fbm_obj$tree
	D    <- fbm_obj$D
	ntips <- fbm_obj$ntips
	
	init_lik = sum(newlnLike(fbm_obj, trait_states, sig2[ind_sig],mu0[ind_mu0],a0[ind_mu0]))
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
				deathRate = calc_deathRates(fbm_obj, trait_states, sig2, ind_sig, mu0, ind_mu0,a0, IND_edge,desc_index_list,1)
		}else if (length(unique(ind_mu0))>1 && death_sig2==2){
				deathRate = calc_deathRates(fbm_obj, trait_states, sig2, ind_sig, mu0, ind_mu0,a0, IND_edge,desc_index_list,2)
		}
		
		deltaRate=sum(deathRate)
		# print(deltaRate)
		
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
						#sig2=c(sig2,sig2[1]) # sample from rate prior
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
						a0_original = a0
						mu0=c(mu0,rnorm(1,0,1)) # sample from rate prior
						a0 =c(a0,rnorm(1,0,0.01))
						# ind_sig in the original order
						ind_mu0=multi_ind_new[new_ind]
						if (length(mu0)>length(unique(ind_mu0))){
							# the new index overlaps completely with the previous
							# remove unused mu0
						 	mu0= mu0_original
							ind_mu0= ind_mu0_original
							a0 = a0_original
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
						a0 = a0[-i]				
						# rescale indexes
						ind_mu0[ind_mu0>i]=ind_mu0[ind_mu0>i]-1					
					}

			}
			
		}
		
	}
	#__ print(list(sig2,ind_sig))	
	return (list(sig2,ind_sig,mu0,a0,ind_mu0))
}

################################## START MCMC ###################################
run_mcmc <- function (fbm_obj,ngen = 100000, control = list(),useVCV=F, sample_f=250,
	                logfile="mcmc.log",update_sig_freq=0.5,dynamicPlot = F,
	                bdmcmc_freq=0.75,useTrend=T,print_freq=100,constRate=F,linTrend=F,
					per_branch_parameters=TRUE, log_anc_states=TRUE, 
                    update_mu0=c(), estimate_HP=FALSE, estimate_HPmu0=FALSE, rate_trend_hp=1,
                    init_mu0=c(), verbose=FALSE){
					
	tree <-      fbm_obj$tree
	x <-         fbm_obj$data
	D <-         fbm_obj$D
	prior_tbl <- fbm_obj$prior_tbl
	ntips <- fbm_obj$ntips 
	PartitionFile <- fbm_obj$PartitionFile
	
	TE=as.matrix(tree$edge)
	a <- mean(fbm_obj$data,na.rm=T)
	y <- rep(a, tree$Nnode - 1)
	sig2 <- c(0.05)  
	mu0 <- c(0)
	a0 <- c(0)
	ind_sig2 <- rep(1, length(c(x, y)))
	ind_mu0  <- rep(1, length(c(x, y)))
	mean_sig <- rep(0, length(c(x, y)))
	mean_anc <- rep(0, length(c(a, y)))
	
	#_ ngen = 2000
	#_ control = list()
	#_ useVCV=F
	#_ sample_f=250
	#_ logfile="mcmc.log"
	#_ update_sig_freq=0.5
	#_ dynamicPlot = F
	#_ bdmcmc_freq=0.8
	#_ useTrend=T
	#_ print_freq=100
	#_ constRate=F
	#_ linTrend=F
	
	if (dynamicPlot == T){dev.new()}
	
	if (PartitionFile != ""){
		bdmcmc_freq = 0
		res_part <- set_model_partitions(fbm_obj,ind_sig2,ind_mu0)
		sig2      <- res_part[[1]]
		mu0       <- res_part[[2]]
		a0        <- res_part[[3]]
		ind_sig2  <- res_part[[4]]
		ind_mu0   <- res_part[[5]]	
		part_ids <- sort(unique(ind_sig2))	
        part_mu_ids <- part_ids
	}


    if (length(fbm_obj$StateTbl) > 0){
        r_tmp <- set_model_trait_partitions(fbm_obj)
        mu0 <- r_tmp[[1]]
        # mu0 = c(0, 0.1, 0.2, 0.3)
        a0 <- r_tmp[[2]]
        ind_mu0 <- r_tmp[[3]]
        part_mu_ids <- sort(unique(ind_mu0))	
        
        # print(head(ind_mu0))
    }
    
    if (length(init_mu0) > 0){
        mu0 = init_mu0
    }
    
	#names(ind_sig2) = c(names(x), D[-1,1])
    if (verbose){
        print(table(ind_sig2))
        print(table(ind_mu0))    
    }
    
    x_tmp = as.vector(table(ind_mu0))
    scaler_w = 1/(1 + log(x_tmp- min(x_tmp) +1))
    std_mu0 <- sd(fbm_obj$data, na.rm=T)*0.025 * scaler_w
    rate_sig2 <- 0.5
    # print(c("std_mu0", std_mu0))
        
    fbm_obj$ind_sig2 = ind_sig2
    fbm_obj$ind_mu0 = ind_mu0
    

	x <- x[tree$tip.label]
	if (is.null(names(y))) {
		names(y) <- length(tree$tip) + 2:tree$Nnode
	}else {y[as.character(length(tree$tip) + 2:tree$Nnode)]}

	# list NA taxa
	ind_NA_taxa = which(is.na(x))
	ind_NAtaxa_in_D               = c()
	brl_NAtaxa_in_D               = c()
	anc_node_of_NAtaxa            = c()
	anc_state_anc_node_of_NAtaxa  = c()
	
	if (per_branch_parameters){
		out_tmp = c("it", "posterior","likelihood","prior", "sig2", "mu0","a0","K_sig2","K_mu0", 
			paste("sig2", 1:length(TE[,1]),sep="_"),paste("mu0", 1:length(TE[,1]),sep="_"),
			paste("a0", 1:length(TE[,1]),sep="_")) 
		}else{
			out_tmp = c("it", "posterior","likelihood","prior", "sig2", "mu0","a0")
			if (PartitionFile != ""){
				out_tmp = c(out_tmp,
				paste("sig2", part_ids, sep="_"),paste("mu0", part_mu_ids,sep="_"),
				paste("a0", part_mu_ids,sep="_"))
			}
		}
	if (log_anc_states){
		out_tmp = c(out_tmp, paste("anc",D[,1],sep="_"), tree$tip.label[ind_NA_taxa])
	}else{
		out_tmp = c(out_tmp, "root", tree$tip.label[ind_NA_taxa])
	}
    if (estimate_HP){
        out_tmp <- c(out_tmp, "sig2_hp")
    }
    if (estimate_HPmu0){
        out_tmp <- c(out_tmp, "mu0_hp")        
    }
    if (logfile != 0){
        cat(c(out_tmp, "\n"), sep="\t", file=logfile, append=F)
    }
	
	
	
	for (tax in 1:length(ind_NA_taxa)){
		ind_NAtaxa_in_D = c(ind_NAtaxa_in_D,   which(D[,2] == ind_NA_taxa[tax]),      which(D[,3] == ind_NA_taxa[tax]))
		brl_NAtaxa_in_D = c(brl_NAtaxa_in_D, D[which(D[,2] == ind_NA_taxa[tax]),4], D[which(D[,3] == ind_NA_taxa[tax]),5])	
		anc_node_of_NAtaxa = c(anc_node_of_NAtaxa, D[ind_NAtaxa_in_D[tax],1])
		anc_state_anc_node_of_NAtaxa = c(anc_state_anc_node_of_NAtaxa, y[which(as.numeric(names(y)) == anc_node_of_NAtaxa[tax])])
	}
	
	# init NAs values
	# phylo imputation
	x[ind_NA_taxa]  = NA
	vector_tip_root_nodes_values = c(x, a, y)
	x_imputed = phylo_imputation(tree, vector_tip_root_nodes_values, sig2[ind_sig2],D,mu0[ind_mu0],a0[ind_mu0],ind_NAtaxa_in_D,
		   				anc_node_of_NAtaxa,brl_NAtaxa_in_D,fbm_obj$dist_from_midpoint,get_expected=1)
	x[ind_NA_taxa]  = x_imputed
	
	fbm_obj$data <- x
	y = runGibbs(fbm_obj,sig2[ind_sig2], c(x, a, y),mu0[ind_mu0],a0[ind_mu0],get_expected=1)[-1]
	
	L   <- newlnLike(fbm_obj, c(x, a, y), sig2[ind_sig2],mu0[ind_mu0],a0[ind_mu0])
    # print(c("LIKELIHOOD:", sum(L), length(L), a, y[1:3], sig2, mu0,sum(ind_sig2), sum(ind_mu0)))
    sd_mu0 = 0.1 
	Pr <- calc_prior(sig2, a, y,mu0, a0, prior_tbl, rate_sig2, sd_mu0)
    if (estimate_HPmu0){
        Pr <- Pr + dexp(sd_mu0, rate_trend_hp, log=T) # add HP probability
    }
    if (estimate_HP){
        Pr <- Pr + dexp(rate_sig2, 1, log=T)
    }
	
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

    if (length(update_mu0) == 0){
        update_mu0 <- 1:length(mu0)
    }

	# START MCMC
	for (i in 1:ngen) {
    
		y.prime     = y
		L.prime     = L
		Pr.prime    = Pr
		a.prime     = a
		sig2.prime  = sig2
		mu0.prime   = mu0
		a0.prime    = a0
        sd_mu0.prime = sd_mu0
        rate_sig2.prime = rate_sig2
		gibbs=0
		hastings=0
    	    	
		# SCREEN OUTPUT / DYNAMIC PLOTS
		if (i %% print_freq==0){
			cat(c("\n",i,round(sum(L),2),length(sig2),length(mu0),a0,a0.prime))
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
	    		        plot.phylo(tree, edge.width=rates_temp[IND_edge]*3, main=paste("lik:",round(sum(L),2),
				  "cat:",length(unique(ind_sig2))),show.tip.label = F)
				
				
			}
		}
    
		j <- (i - 1)%%(tree$Nnode + 1)
		rr= runif(2,0,1)
		
		update_sig_freq=0.95
		if (rr[1]<update_sig_freq) {
			if (rr[2]< 0.33 && i > 50){ # SIG2 UPDATE
				s_ind  = sample(1:length(sig2),1)
				sig2_update <-  update_multiplier_proposal(sig2.prime[s_ind],1.05)
				sig2.prime[s_ind] = sig2_update[1]
				hastings = sig2_update[2]			
			}
			else if (rr[2]<0.66 && useTrend==T){ # MU0 UPDATE
                if (length(update_mu0) == 1){
                    m_ind = update_mu0[1]
                }else{
                    m_ind  = sample(update_mu0,1) #<---- ONLY UPDATE thise within update_mu0
                }
				mu0_update <-  mu0[m_ind] + rnorm(n = 1, sd = std_mu0[m_ind])
				if (linTrend){
					# update diff
					# r_diff = mu0[m_ind] - a0[m_ind]
					# r_diff_update <- r_diff + rnorm(n = 1, sd = sd(fbm_obj$data, na.rm=T)*0.0005)
					# a0_update <- mu0_update - r_diff_update
					# print(c(a0[m_ind], a0_update, r_diff_update, r_diff))
					a0_update <-  a0[m_ind] + rnorm(n = 1, sd = sd(fbm_obj$data, na.rm=T)*0.0005)
					a0.prime[m_ind] = a0_update
				}	
				mu0.prime[m_ind] = mu0_update
				
			}
			else{ # ROOT STATE
				a.prime <- a + rnorm(n = 1, sd = 0.5) #sqrt(con$prop[j + 1]))
                if (estimate_HP){
                    rate_sig2_update <- update_multiplier_proposal(rate_sig2, 1.05)
                    rate_sig2.prime <- rate_sig2_update[1]
                    hastings <- hastings + rate_sig2_update[2]
                }
                if (estimate_HPmu0){
    				sd_mu0_update <-  update_multiplier_proposal(sd_mu0, 1.05)
    				sd_mu0.prime <- sd_mu0_update[1]
                    hastings <- hastings + sd_mu0_update[2]
                }
			}
			
		}
		 else{
			# ANC STATES
			#useTrend=F
			RUNIF = runif(3,0,1)
			if (RUNIF[1]<bdmcmc_freq && i > 100){
				trait_states=c(x, a.prime, y.prime)
				if (useTrend){
					death_sig2=sample(2,1)
				}else{death_sig2=1}
				if (constRate==T){death_sig2=2}	
					bdmcmc_list   = run_BDMCMC(fbm_obj, trait_states, sig2.prime,ind_sig2, mu0.prime, ind_mu0, a0.prime, 
										IND_edge,new_ind,desc_index_list,death_sig2)
					sig2.prime    = bdmcmc_list[[1]]
					ind_sig2      = bdmcmc_list[[2]]					
					mu0.prime     = bdmcmc_list[[3]]
					a0.prime      = bdmcmc_list[[4]]
					ind_mu0       = bdmcmc_list[[5]]										
					#print (ind_mu0)
					gibbs=1
			} else {
				
				if (runif(1)<0.1){
					# phylo imputation
					x[ind_NA_taxa]  = NA
					vector_tip_root_nodes_values = c(x, a.prime, y.prime)
					x_imputed = phylo_imputation(tree, vector_tip_root_nodes_values, sig2.prime[ind_sig2],D,mu0.prime[ind_mu0],a0.prime[ind_mu0],
										ind_NAtaxa_in_D,anc_node_of_NAtaxa,brl_NAtaxa_in_D,fbm_obj$dist_from_midpoint,get_expected=0)
					x[ind_NA_taxa]  = x_imputed
					fbm_obj$data <- x
					gibbs=1
				}else{
					# gibbs update anc states
					vector_tip_root_nodes_values = c(x, a.prime, y.prime)
					y_temp = runGibbs(fbm_obj,sig2.prime[ind_sig2], vector_tip_root_nodes_values,mu0.prime[ind_mu0],a0.prime[ind_mu0],get_expected=0)
					a.prime= y_temp[1]
					y.prime= y_temp[-1]
					#x[ind_NA_taxa]  = NA
					#vector_tip_root_nodes_values = c(x, a.prime, y.prime)
					#x_imputed = phylo_imputation(tree, vector_tip_root_nodes_values, sig2.prime[ind_sig2],D,mu0.prime[ind_mu0],a0.prime[ind_mu0],
					#					ind_NAtaxa_in_D,anc_node_of_NAtaxa,brl_NAtaxa_in_D,fbm_obj$dist_from_midpoint,get_expected=0)
					#x[ind_NA_taxa]  = x_imputed
					#fbm_obj$data <- x
					
					gibbs=1					
				}
				
			}
		}
		
		
		
		
	       # calc post
		L.prime <- newlnLike(fbm_obj, c(x, a.prime, y.prime), sig2.prime[ind_sig2],mu0.prime[ind_mu0],a0.prime[ind_mu0])			
		Pr.prime <- calc_prior(sig2.prime, a.prime, y.prime, mu0.prime, a0.prime, prior_tbl, rate_sig2.prime, sd_mu0.prime)
        if (estimate_HPmu0){
            Pr.prime <- Pr.prime + dexp(sd_mu0.prime,
                                        rate=rate_trend_hp, log=T)
        }
        if (estimate_HP){
            Pr.prime <- Pr.prime + dexp(rate_sig2.prime, 1, log=T)            
        }
		
        # print( c(sum(L.prime), sum(Pr.prime), sig2.prime, a.prime, mu0.prime, x[ind_NA_taxa], sum(Pr), sum(L)) )
		
		
		if ( (sum(Pr.prime) + sum(L.prime) - sum(Pr) - sum(L) + hastings) >= log(runif(1,0,1)) || gibbs==1){    
			y    = y.prime
			L    = L.prime
			Pr   = Pr.prime
			a    = a.prime
			sig2 = sig2.prime
			mu0  = mu0.prime
			a0   = a0.prime
            sd_mu0 = sd_mu0.prime
            rate_sig2 = rate_sig2.prime
		}
 
	     if (i %% sample_f == 0) {
 			rates_temp=sig2[ind_sig2]
 			trends_temp=mu0[ind_mu0]
 			trend_trends = a0[ind_mu0]
			
			if (log_anc_states){
				anc_tmp = c(a, y)
			}else{
				anc_tmp = c(a)
			}
            x_imputed_looged = x_imputed
            if (estimate_HP){
                x_imputed_looged = c(x_imputed_looged, rate_sig2)
            }
            if (estimate_HPmu0){
                x_imputed_looged = c(x_imputed_looged, sd_mu0)
            }
            
			
			 if (per_branch_parameters){
                 if (logfile != 0){
    	 			cat(c(i,sum(L)+sum(Pr), sum(L),sum(Pr), mean(rates_temp),mean(trends_temp),mean(trend_trends), length(sig2),length(mu0), 
    	 				rates_temp[IND_edge],trends_temp[IND_edge],trend_trends[IND_edge], anc_tmp, x_imputed_looged, "\n"),sep="\t", file=logfile, append=T)
                    }
			 }else{
				if (PartitionFile != ""){
					prm_tmp = c(sig2, mu0, a0)
				}else{
					prm_tmp = c()
				}
                if (logfile != 0){
     	 			cat(c(i,sum(L)+sum(Pr), sum(L),sum(Pr), mean(rates_temp),mean(trends_temp),mean(trend_trends), prm_tmp,
     	 				anc_tmp, x_imputed_looged, "\n"),sep="\t", file=logfile, append=T) 			 	                    
                }
			 }
			#print(x[ind_NA_taxa])
	    }
    
	}
    return (fbm_obj)
}

################################## END   MCMC ###################################


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



plot_results <- function(fbm_obj, logfile, resfile="results.pdf" , exp_trait_data=0, return_estimated_prm=FALSE, width=15, height=15){	
	require(phytools)
	require(methods)
	library(scales)
	library(plotrix)
	out_tbl = read.table(logfile,header=T)
	ind_sig2_col = grep('sig2_', colnames(out_tbl), value=F)
	burnin= round(0.25*dim(out_tbl)[1])

	pdf(file=resfile,width=width, height=height)

	tree <- fbm_obj$tree
	ntips <- fbm_obj$ntips
	data <- fbm_obj$data
	rescale_trait_data <- fbm_obj$trait_rescaling
	
	sp_with_missing_data = names(data[which(is.na(data))])
	phylo_imputation_data <- out_tbl[sp_with_missing_data]
	data[which(is.na(data))] = apply(phylo_imputation_data,FUN=mean,2)
	
	# Rates section

	rates_temp_median = apply(out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=median,2)
	rates_temp_mean = apply(out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=mean,2)
	estK_all =out_tbl$K_sig2[burnin:dim(out_tbl)[1]]
	h_rate=hist(estK_all, breaks=0:max(estK_all)+0.5,plot=F)
	best_K = which(h_rate$counts==max(h_rate$counts))
	estK = round(quantile(estK_all,probs = c(0.025,0.975)))
	rates_min = apply(X=out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=quantile, probs=0.025,2)
	rates_max = apply(X=out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=quantile, probs=0.975,2)

	par(mfcol=c(2,1))

	temp<-log(rates_temp_mean) + max(abs(log(rates_temp_mean))) + 0.5

	# Rates are log-transformed to avoid extremes, then make sure that values are positive, then 
	# add an arbitrary small number to avoid zeros for the cex - edge width parameter in plotting.
	
	col = rep("black",length(tree$edge.length))
	names(col)=c(1:length(tree$edge.length))
	plot.phylo(tree, edge.width=temp, main=paste("Estimated rates,",sprintf("K: %s (%s-%s)",best_K,estK[1],estK[2])),show.tip.label = T, cex=0.1, align.tip.label=T)



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
	if (xa==0){xa=0.01}
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
	plot.phylo(tree, edge.width=1, main=paste("Estimated trends,",sprintf("K: %s (%s-%s)",best_K,estK[1],estK[2])),show.tip.label = T, cex=0.1, edge.color=edge_cols)
	testcol<-rbPal(n=length(breaksList))

	col.labels<-c(format(min_col, digits = 4), "","0","", format(max_col, digits = 4))
  	limit_plot<-max(axisPhylo())
  	color.legend(0,0.15,0+(limit_plot)/10,25,col.labels,rect.col=testcol, gradient="y", cex=0.8, align="lt")

    par(mfrow=c(3,2))
	C = h_rate$counts
	names(C) = 1:max(estK_all)
	barplot(C,main="Estimated number of shifts - Sig2")

	plot(rates_temp_mean,  main=paste("sig2s, ", sprintf("K: %s (%s-%s)",best_K,estK[1],estK[2])),pch=21,col="darkblue",bg="blue", ylim=c(min(rates_min)*0.5,max(rates_max)*1.5),type="n")
	points(rates_temp_mean,pch=21,col="darkred",bg="red")
	segments(x0= 1:length(rates_min), y0=rates_min, x1=1:length(rates_min), y1=rates_max,col=alpha("darkred",0.4),lwd=2)

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
	sampled_indx = sample(1:dim(use_mcmc)[1], size=sample_size, replace = FALSE)
	
	sample_mcmc <- use_mcmc [ sampled_indx , ]
	colnames(sample_mcmc) = seq(from=ntips+1, to=max(tree$edge[,2]), by=1)
	sample_mcmc <- as.matrix(sample_mcmc)
	sample_phylo_imputation_data <- phylo_imputation_data[sampled_indx,]
	

	# Coming back to the original values of the trait, and applying the same transform/scaling used for MCMC
	if(exp_trait_data == 0) {
	   sample_mcmc_raw <- sample_mcmc
	   est_anc_states_mean <- est_anc_states
	   data_raw <- data
	   sample_phylo_imputation_data_raw <- sample_phylo_imputation_data
	  	}
	if(exp_trait_data == 1){ 
	  sample_mcmc_raw <- exp(sample_mcmc)	 
	  est_anc_states_mean <- exp(est_anc_states)
	  data_raw <- exp(data)	
	  sample_phylo_imputation_data_raw <- exp(sample_phylo_imputation_data)
	   }
	if(exp_trait_data == 10){
	  sample_mcmc_raw <- 10^sample_mcmc
	  est_anc_states_mean <- 10^est_anc_states
	  data_raw <- 10^data	
	  sample_phylo_imputation_data_raw <- 10^sample_phylo_imputation_data
    	}		
	
	
	
	
	
	sample_mcmc_raw = sample_mcmc_raw/rescale_trait_data
	est_anc_states_mean_raw <- est_anc_states_mean/rescale_trait_data
 	data_raw <- data_raw/rescale_trait_data
	sample_phylo_imputation_data_raw <- sample_phylo_imputation_data_raw/rescale_trait_data
	# mean values phenogram 

	color=adjustcolor( "darkred", alpha.f = 1)
	est_anc_states_mean_raw[is.na(est_anc_states_mean_raw)] = mean(est_anc_states_mean_raw, na.rm=T)
	phenogram_invTime(tree, est_anc_states_mean_raw, col=color, main="Traigram", spread.labels=F)
 
	# mcmc sample of phenograms

	alphacol = 1/sample_size + sample_size*0.001

	color=adjustcolor( "darkred", alpha.f = alphacol)
	temp = c(data_raw, sample_mcmc_raw[1,])	
	temp[is.na(temp)] = mean(temp, na.rm=T)
	phenogram_invTime(tree, temp, col=color, main="Traigram", spread.labels=F, ftype="off", 
			ylim=c(min(c(sample_mcmc_raw,temp)), max(c(sample_mcmc_raw,temp))))
	
	for (s in 2:dim(sample_mcmc)[1]) { 
		data_tmp <- data_raw
		data_tmp[which(is.na(fbm_obj$data))] <- as.numeric(sample_phylo_imputation_data_raw[sampled_indx[s], ])
		
		temp = c(data_tmp, sample_mcmc_raw[s,])	
		phenogram_invTime(tree, temp, col=color, main="Traigram", spread.labels=F, ftype="off", add=T)
	}

	n<-dev.off()
	if (return_estimated_prm==TRUE){
		
		ind_a0_col = grep('a0_', colnames(out_tbl), value=F)
		a0_temp_mean = (apply(out_tbl[burnin:dim(out_tbl)[1], ind_a0_col],FUN=mean,2))
		a0_min = apply(X=out_tbl[burnin:dim(out_tbl)[1], ind_a0_col],FUN=quantile, probs=0.025,2)
		a0_max = apply(X=out_tbl[burnin:dim(out_tbl)[1], ind_a0_col],FUN=quantile, probs=0.975,2)
		
		ind_anc_col = grep('anc_', colnames(out_tbl), value=F)[1]
		root_state_mean = mean(out_tbl[burnin:dim(out_tbl)[1], ind_anc_col])
		root_state_min = quantile(out_tbl[burnin:dim(out_tbl)[1], ind_anc_col], probs=0.025,2)
		root_state_max = quantile(out_tbl[burnin:dim(out_tbl)[1], ind_anc_col], probs=0.975,2)
		
		
		prm = c(unique(rates_temp_mean), unique(rates_min),unique(rates_max),
		unique(mu0_temp_mean), unique(mu0_min), unique(mu0_max), 
		unique(a0_temp_mean), unique(a0_min), unique(a0_max),
		root_state_mean,root_state_min,root_state_max)
		
		res = list()
		res[[1]] <- est_anc_states_mean_raw
		res[[2]] <- prm
		
		return(res)
	}
		

}




plot_res_trend <-function(res,time_axis, ylab,main=""){
	time_axis = sort(-time_axis)
	library(HDInterval)
	CI = hdi(res)
	CI_max = CI["upper",]
	CI_min = CI["lower",]
	significantly_positive = time_axis[which(CI_min > 0)]
	significantly_negative = time_axis[which(CI_max < 0)]
	plot(time_axis,res[1,], ylab=ylab,xlab="Time", type="n",
		 ylim = c(min(res), max(res)), main=main)
	abline(0,0,lty=2)
	points(significantly_positive, rep(min(res), length(significantly_positive)), col="darkred")
	points(significantly_negative, rep(min(res), length(significantly_negative)), col="blue")
	
	for (s in 1:min(100,dim(res)[1])){
		lines(time_axis,res[s,],col=alpha("black",0.2))
	}
	CIm = t(as.data.frame(CI))
	return(CIm)
}




plot_time_varying_trend <- function(fbm_obj, logfile, resfile="trends.pdf", height=15, width=10){
	out_tbl = read.table(logfile,header=T)
	ind_mu0_col = grep('mu0_', colnames(out_tbl), value=F)
	ind_a0_col = grep('a0_', colnames(out_tbl), value=F)
	burnin= round(0.25*dim(out_tbl)[1])
	out_tbl <- out_tbl[burnin:dim(out_tbl)[1], ]
	
	root_age = max(fbm_obj$dist_from_the_root)
	
	mu0_temp_mean = unique(apply(out_tbl[, ind_mu0_col],FUN=mean,2))
	a0_temp_mean = unique(apply(out_tbl[, ind_a0_col],FUN=mean,2))
	
	pdf(file=resfile,width=width, height=height)
	par(mfrow = c(length(mu0_temp_mean),2))
	
	t = seq(-root_age*0.5,root_age*0.5,length.out=1000)
	time_axis = seq(0,root_age,length.out=1000)
	for (i in 1:length(mu0_temp_mean)){
		
		res = NULL
		for (s in 1:1000){
			indx = sample(1:dim(out_tbl)[1], size=1)
			mu0_s = unique(as.numeric(out_tbl[indx, ind_mu0_col]))
			a0_s = unique(as.numeric(out_tbl[indx, ind_a0_col]))
			evolving_trend = a0_s[i]*t*t + mu0_s[i]*t + 0
			evolving_trend = (evolving_trend - evolving_trend[length(evolving_trend)]) *(1/fbm_obj$trait_rescaling)
			res = rbind(res,evolving_trend)
		}
		sig_res_phenotype <- plot_res_trend(res, time_axis, ylab="Change in expected phenotype (y_t)", 
					main = paste0("Partition ", i))
		mean_res = apply(res,FUN=mean,2)
		lines(sort(-time_axis),mean_res, lwd=2,col="red")
		if (i == 1){
			res_tbl = cbind(-time_axis, mean_res, sig_res_phenotype)
			colnames(res_tbl) <- c("time", paste0("Partition-", i, "_mean"), paste0("Partition-", i, "_min"), paste0("Partition_", i, "-max"))
		} else{
			temp = cbind(mean_res, sig_res_phenotype)
			colnames(temp) <- c(paste0("Partition_", i, "_mean"), paste0("Partition_", i, "_min"), paste0("Partition_", i, "_max"))
			res_tbl = cbind(res_tbl, temp)
		}
		
		
		
		res = NULL
		for (s in 1:1000){
			indx = sample(1:dim(out_tbl)[1], size=1)
			mu0_s = unique(as.numeric(out_tbl[indx, ind_mu0_col]))
			a0_s = unique(as.numeric(out_tbl[indx, ind_a0_col]))
			evolution_of_the_trend = a0_s[i]*t + mu0_s[i]
			res = rbind(res,evolution_of_the_trend)
		}
		sig_res_trend <- plot_res_trend(res, time_axis, ylab="Trend parameter (mu_t)")
		mean_res = apply(res,FUN=mean,2)
		lines(sort(-time_axis),mean_res, lwd=2,col="red")
		
		if (i == 1){
			res_trend = cbind(-time_axis, mean_res, sig_res_phenotype)
			colnames(res_tbl) <- c("time", paste0("Partition-", i, "_mean"), paste0("Partition-", i, "_min"), paste0("Partition_", i, "-max"))
		} else{
			temp = cbind(mean_res, sig_res_phenotype)
			colnames(temp) <- c(paste0("Partition_", i, "_mean"), paste0("Partition_", i, "_min"), paste0("Partition_", i, "_max"))
			res_trend = cbind(res_trend, temp)
		}
		
		
		
	}
	n<- dev.off()
	
	return(list(res_tbl, res_trend))
	
}


get_r2_mse <-function(x,y){
	r2 = summary(lm(y~x))$r.squared
	mse = mean((y-x)^2) # mean square error
	return(c(r2, mse))
}



simulate_trait_data <- function(fbm_obj, sigma2=0.2, mu0=0, a0=0, 
                                plot=F, seed=0, plot_file="",color_by_mu0=T,
                                randomize_mu0=NULL
                                ){
	if (seed > 0){
		set.seed(seed)
	}
	ntips	  <- fbm_obj$ntips
	tree <- fbm_obj$tree
	D		  <- fbm_obj$D
	prior_tbl  <- fbm_obj$prior_tbl
	root_dist <- fbm_obj$dist_from_the_root
    dist_from_midpoint <- fbm_obj$dist_from_the_root
	
	root_state = 0
	all_states = rep(root_state, (ntips*2-1))
	root_age = -max(root_dist)
    
    # rate shifts
    sigma2 <- sigma2[fbm_obj$ind_sig2]
    mu0 <- mu0[fbm_obj$ind_mu0]
    a0 <- a0[fbm_obj$ind_mu0]
    
    if (length(randomize_mu0) > 0){
        for (i in 1:length(unique(fbm_obj$ind_mu0))){
            ind_tmp <- which(fbm_obj$ind_mu0 == i)
            mu0[ind_tmp] <- mu0[ind_tmp] + rnorm(length(ind_tmp), 0, randomize_mu0[i])
        }
    }
    
	
	for (indx in (ntips+1):(ntips*2-1)){
		
		i = which(D[,1]==indx)
		
		anc_ind <- D[i,1];
		a_ind   <- D[i,2];  # index of descendants
		b_ind   <- D[i,3];  # index of descendants
		vpa	 <- D[i,4];  # br length
		vpb	 <- D[i,5];  # br length
		anc = all_states[anc_ind]
        if (is.na(sum(all_states))){break}
    	anc_ind_v <- D[i,1]
    	a_ind_v   <- D[i,2]
    	b_ind_v   <- D[i,3]
        
        if (a_ind_v > ntips){
            a_ind_v <- a_ind_v -1
        }
        if (b_ind_v > ntips){
            b_ind_v <- b_ind_v -1
        }
        
        
		m1 <- anc + (vpa*mu0[a_ind_v] + a0[a_ind_v]*vpa*dist_from_midpoint[a_ind])
		s1 <- sqrt((vpa)*sigma2[a_ind_v])
		m2 <- anc + (vpb*mu0[b_ind_v] + a0[b_ind_v]*vpb*dist_from_midpoint[b_ind])
		s2 <- sqrt((vpb)*sigma2[b_ind_v])
		
        # cat(c("\n", round(indx), round(i),  m1, s1, m2, s2, "\n"))
        # print(all_states)
		all_states[a_ind] <- rnorm(1, m1, s1)
		all_states[b_ind] <- rnorm(1, m2, s2)
				
	}
    
	if (plot){
		if (plot_file != ""){
			pdf(file=plot_file, 18,18)
		}
        if (sum(a0) > 0){
    		par(mar=c(2,4,1,1)*1.2)
    		par(fig=c(0,10,4,10)/10)            
        }
        
        if (color_by_mu0){
        	rbPal <- colorRampPalette(c("#2166ac","gray","#b2182b"))
        	# maximum color scale for trend 
        	xa <- max(abs(c (mu0)))+ (0.1* max(abs(c (mu0))))
        	if (xa==0){xa=0.01}
        	min_col <- -xa
        	max_col <- xa
        	beta_shape <- 0.5
        	test<- qbeta(seq(0,1, by=0.05), beta_shape, beta_shape, lower.tail = TRUE, log.p = FALSE)
        	x = test * (max_col-min_col) 
        	breaksList = x + min_col
        	edge_cols <- rbPal(length(breaksList))[as.numeric(cut(sort(mu0),breaks = breaksList))]
        	names(edge_cols) <- order(mu0)
        	edge_cols <- edge_cols [order(as.numeric(names(edge_cols)))]
        	names(edge_cols)<- NULL
            
        }
        
        
		# plot(-(max(root_dist)-root_dist)[1:fbm_obj$ntips], all_states[1:fbm_obj$ntips],)
		plot(NULL, xlim=c(root_age, 0), ylim=c(min(all_states), max(all_states)),
			xlab="Time", ylab="Trait")
         
		for (indx in (ntips+1):(ntips*2-1)){	
			i = which(D[,1]==indx)
			anc_ind <- D[i,1];
			a_ind   <- D[i,2];  # index of descendants
			b_ind   <- D[i,3];  # index of descendants
			vpa	 <- D[i,4];  # br length
			vpb	 <- D[i,5];  # br length
            
			anc = all_states[anc_ind]
			x0 = root_age + root_dist[anc_ind]
			x1 = root_age + root_dist[a_ind]
			segments(x0, anc, x1, all_states[a_ind], col=edge_cols[a_ind])
			x1 = root_age + root_dist[b_ind]
			segments(x0, anc, x1, all_states[b_ind], col=edge_cols[b_ind])
		
		}
        if (sum(a0) > 0){
    		t = seq(-root_age*0.5,root_age*0.5,length.out=100) - (root_age/2)
    		time_axis = seq(0,root_age,length.out=100)
		
    		trend = a0*t + mu0
    		par(mar=c(4,4,1,1)*1.2)
    		par(fig=c(0,10,0,4)/10)
    		par(new=T)
		
    		plot(NULL, xlim=c(root_age, 0), ylim=c(-1, 1), xlab="Time", ylab="Trend")
    		points(time_axis, trend, type="l")
    		abline(h=0,lty=2, col="red")
		}
		if (plot_file != ""){
			dev.off()
		}
		
	}
    tip_states = all_states[1:ntips]
    names(tip_states) <- names(fbm_obj$data)
    res <- NULL
    res$tip_states <- tip_states
    res$all_states <- all_states
    res$sigma2 <- sigma2
    res$mu0 <- mu0
    res$a0 <- a0
	return(res)
}	
	