#!/usr/bin/Rscript
# version 20150713

arg <- commandArgs(trailingOnly=TRUE)

require(phytools)
require(diversitree)
library(gtools)

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
	#return((-(xa - xb)^2 / (2 * sigma2 * (vpa + vpb)) - log(sqrt(2 * pi * sigma2 * (vpa + vpb)))));
}

fN2<-function(xa, xb, vpa, vpb, sigma2_a,sigma2_b, anc) {
	#same as fN but with a known ancestor instead of xa-xb, see Joe's book eqn 23.10 (some mistakes in the denominator though in his equation)
	return(dnorm(xa, mean=anc, sd= sqrt(vpa*sigma2_a),log=T) + dnorm(xb, mean=anc, sd= sqrt(vpb*sigma2_b),log=T)  )
	#return(-(((xa - anc)^2 / vpa) + ((xb - anc)^2 / vpb)) / (2 * sigma2) - log(sqrt(2 * pi * sigma2 * (vpa + vpb))));
}

newlnLike<-function(tree, vector_tip_root_nodes_values, sigma2,D) {
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
	
        L1 = dnorm(a_v, mean=anc_v, sd= sqrt(vpa_v*sigma2[a_ind_v]),log=T)
        L2 = dnorm(b_v, mean=anc_v, sd= sqrt(vpb_v*sigma2[b_ind_v]),log=T)  
         
        L = L1+L2
	return(L) 
	
	#print(condvec)
	#return(lik) 
  
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

runGibbs <- function(sigma2, vector_tip_root_nodes_values,D,prior_tbl,get_expected=0) {
	
	#prior_tbl = calc_prior_dist_root_calibration(D,sigma2)
	
	vec_values = vector_tip_root_nodes_values
        # loop over Johnatan's tbl from most recent to root
        for (i in dim(D)[1]:2){
	#for (rep in 1:20){
	#	i = sample(dim(D)[1])[1]
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

		

		desc_lik = get_joint_mu_s2(a, (vpa)*sigma2[sig_a_ind], b, (vpb)*sigma2[sig_b_ind])
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
			
			#if (stem_ind>ntips){
			sig_stem_ind=stem_ind-1
			#}else{sig_stem_ind=stem_ind}
			
			
			stem_val = vec_values[stem_ind]		
			stem_lik = c(stem_val, (stem_brl)*sigma2[sig_stem_ind])
			
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


################ DPP GIBBS SAMPLER
G0 <- function (){
	return(abs(rcauchy(1,0,1)))
}

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

cond_alpha_proposal <- function(current_alpha,k,n,hp_gamma_shape=1,hp_gamma_rate=0.2){
	z = c(current_alpha + 1.0, n)
	f = rdirichlet(1, z)
	eta = f[1]	
	u = runif(1) 
	x = (hp_gamma_shape + k - 1.0) / ((hp_gamma_rate - log(eta)) * n)
	if ((u / (1.0-u)) < x){
		new_alpha = rgamma(1, shape=(hp_gamma_shape+k), scale=(1./(hp_gamma_rate-log(eta))) )
	}
	else{
		new_alpha = rgamma(1, shape=(hp_gamma_shape+k-1.), scale=1./(hp_gamma_rate-log(eta)) )
	}
	return (new_alpha)	
}

DDP_gibbs_sampler4 <- function(tree,sigma2,ind,x,a,y,D,prior_tbl,alpha_par_Dir,ind_anc,des_col){
	vec_values=c(x,a,y)
	ay_temp = c(a,y)
	# sig2: sig2ameters for each category
	n_data=length(ind)
	# GIBBS SAMPLER for NUMBER OF CATEGORIES - Algorithm 4. (Neal 2000)
	sig2=sigma2 # sig2   parameters for each category

	eta=rep(0,length(sig2)) # number of elements in each category
	for (j in 1:length(sig2)){
		eta[j] = length(ind[ind==j])
	}
	#__ cat("\nstart ind",ind)
	#__ cat("\nstart eta",eta)
	
	u1 = runif(n_data,0,1) # init random numbers
	new_lik_vec=rep(0,n_data) # store new sampled likelihoods
	
	new_alpha_par_Dir = cond_alpha_proposal(alpha_par_Dir,length(sig2),n_data,hp_gamma_shape=1,hp_gamma_rate=0.2)
	#__ cat("\n\nLIK1:", sum(newlnLike(tree, vec_values, sig2[ind],D)))
	
	i = ind_anc
	
	anc_ind <- D[i,1];
	a_ind   <- D[i,2];  # index of descendants
	b_ind   <- D[i,3];  # index of descendants
	vpa     <- D[i,4];  # br length
	vpb     <- D[i,5];  # br length
	anc = vec_values[anc_ind]
	a = vec_values[a_ind]
	b = vec_values[b_ind]
	
	
	if (a_ind>ntips){
		des_a_ind=a_ind-1 # because no sig2 value for root
	}else{des_a_ind=a_ind}

	if (b_ind>ntips){
		des_b_ind=b_ind-1
	}else{des_b_ind=b_ind}
	

	#__ cat("\nETA",eta)
	#anc_ind <- TE[i,1];
	#des_ind <- TE[i,2];  # index of descendant
	#anc = vec_values[anc_ind]
	#des = vec_values[des_ind]
	
	
	if (des_col==1){
		des_ind= des_a_ind
	}else{des_ind=  des_b_ind}
	
	
	k1 = length(sig2)

	if (length(ind[ind==ind[des_ind]])==1){ # is singleton
		k1 = k1 - 1
		sig2_k1 = sig2			
		if (u1[i]<= k1/(k1+1)){
			0 #NULL
		}else{
			ind[des_ind] = k1 + 1
		} # this way n_ic for singleton is not 0
	} else{ # is not singleton
		sig2_k1 = append(sig2,G0())
	}
	
	# construct prob vector
	lik_vec=c()
	ay=list()
	for (R in sig2_k1){
		sigma2_vec=sig2[ind]
		sigma2_vec[des_ind]=R
		l_temp = sum(newlnLike(tree, vec_values, sigma2_vec,D))
		lik_vec = append(lik_vec, l_temp)
	}

	rel_lik = calc_rel_prob(lik_vec)

	if (length(sig2_k1)>length(eta)){ # sig2_k1 add one element only when i is not singleton
		eta[ind[des_ind]] = eta[ind[des_ind]]-1
		eta_temp= append(eta,new_alpha_par_Dir/(k1+1))
	}
	else{
		eta_temp = eta
	}
	
	P=eta_temp*rel_lik
	IND = random_choice_P(P)[2] 
	ind[des_ind] = IND # update ind vector
	if (IND==(length(sig2_k1))){
		sig2 = sig2_k1 # add category
	}
	# Change the state to contain only those sig2 are now associated with an observation
	# create vector of number of elements per category
	eta=rep(0,length(sig2))
	for (j in 1:length(sig2)){
		eta[j] = length(ind[ind==j])
	}
	#__ cat( c("\n\nnew eta:",eta, "\nsig:",sig2_k1))
	# remove parameters for which there are no elements
	sig2 = sig2[eta>0]
	# rescale indexes
	ind_rm = which(eta==0) # which category has no elements
	if (length(ind_rm)>0){
		ind[ind>=ind_rm] = ind[ind>=ind_rm]-1
		# update eta
		eta = eta[-ind_rm]
	}
	new_lik_vec[des_ind]=lik_vec[IND]

	likA = sum(new_lik_vec) # newlnLike(tree, vector_tip_root_nodes_values, sig2[ind],D)
	sig2A = sig2
	return (list(sig2A, ind,new_alpha_par_Dir))
	
}

################################## START MCMC ###################################
mcmc.gibbs4 <- function (tree, x,alter_ind, D, prior_tbl,true_anc= NA,ngen = 100000, control = list(), gibbs_sampler=T,useVCV=F, sample=100,logfile="log",update_sig_freq=0.5,useDPP=T,dynamicPlot = F) 
{
	TE=as.matrix(tree$edge)
	cat(c("it", "posterior","likelihood","prior", "sig2","K", "root","alpha", paste("sig2", 1:length(TE[,1]),sep="_"), D[-1,1],"\n"),sep="\t", file=logfile, append=F)
	a <- 0
	y <- rep(0, tree$Nnode - 1)
	sig2 <- c(0.2) # rep(0.2, length(c(x,  y))) # 
	ind_sig2 <- rep(1, length(c(x, y))) #1:length(c(x, y)) #
	mean_sig <- rep(0, length(c(x, y)))
	if (dynamicPlot == T){dev.new()}
	#new_alter_ind = tree$edge[alter_ind,2]
	#new_alter_ind[new_alter_ind>ntips] = new_alter_ind[new_alter_ind>ntips]-1
	#ind_sig2[new_alter_ind] = 2
	#names(ind_sig2) = c(names(x), D[-1,1])
	alpha_par_Dir = 1

	# lik function
	if (useVCV ==T){
		temp <- phyl.vcv(as.matrix(x), vcv(tree), 1)
		likelihood <- function(C, invC, detC, x, sig2, a, y) {
			z <- c(x, y) - a
			logLik <- -z %*% invC %*% z/(2 * sig2) - nrow(C) * log(2 * pi)/2 - nrow(C) * log(sig2)/2 - detC/2
			return(logLik)
		}
		C <- vcvPhylo(tree)
		if (any(tree$edge.length <= (10 * .Machine$double.eps))) 
		stop("some branch lengths are 0 or nearly zero")
		invC <- solve(C)
		detC <- determinant(C, logarithm = TRUE)$modulus[1]    	
	}
	else{
	    likelihood <- function(C, invC, detC, x, sig2, a, y) {
		vector_tip_root_nodes_values = c(x, a, y)
		return(sum(newlnLike(tree, vector_tip_root_nodes_values, sig2,D)))
		}    	
	}

	# priors
	log.prior <- function(sig2, a, y) {
		prior_sig2 <- sum(dcauchy(sig2, loc=0,scale=1, log = TRUE) )
		# + sum(dnorm(c(a, y), mean = pr.mean[2:length(pr.mean)], sd = sqrt(pr.var[1 + 1:tree$Nnode]), log = TRUE))
		# prior_tbl = calc_prior_dist_root_calibration(D,sig2)
		prior_root = sum(dnorm(c(a), mean = prior_tbl[,1], sd = prior_tbl[,2], log = T))
		prior_anc = sum(dnorm(c(a, y), mean = 0, sd = 100, log = T))
		return(prior_sig2+prior_root+prior_anc)
	}

	x <- x[tree$tip.label]
	if (is.null(names(y))) {
		names(y) <- length(tree$tip) + 2:tree$Nnode
	}
	else {y[as.character(length(tree$tip) + 2:tree$Nnode)]}

	L  <- newlnLike(tree, c(x, a, y), sig2[ind_sig2],D)
	Pr <- log.prior(sig2, a, y)


	IND_edge = c()
        for (ind_edge in tree$edge[,2]){
		if (ind_edge>ntips){
			ind_edge_t=ind_edge-1
		}else{ind_edge_t=ind_edge}
        	IND_edge = append(IND_edge,ind_edge_t)
        }



	# START MCMC
	for (i in 1:ngen) {
    
		y.prime     = y
		L.prime     = L
		Pr.prime    = Pr
		a.prime     = a
		sig2.prime  = sig2
		gibbs=0
		hastings=0
    	    	
		
		if (dynamicPlot == T){
			mean_sig = rbind(mean_sig,sig2[ind_sig2])
		
			if (i%%print_f == 0 || i==1) {
				cat(c("\n",i,round(c(sum(L),sig2),2),"\n",ind_sig2))
				start = round(dim(mean_sig)[1]*0.1)
				end = dim(mean_sig)[1]
	    		        rates_temp = apply(mean_sig[start:end,], 2, FUN=mean)
	    		        plot.phylo(tree, edge.width=rates_temp[IND_edge]*3, main=paste("lik:",round(sum(L),2),"cat:",length(sig2)),show.tip.label = F)
			}
		}
    
		j <- (i - 1)%%(tree$Nnode + 1)
		rr= runif(1,0,1)
		
		update_sig_freq=0.9
		if (rr<update_sig_freq) {
			if (runif(1,0,1)>0.5){
				s_ind  = sample(1:length(sig2),1)
				sig2_update <-  update_multiplier_proposal(sig2.prime[s_ind],1.2)
				sig2.prime[s_ind] = sig2_update[1]
				hastings = sig2_update[2]			
			}
			else{
				a.prime <- a + rnorm(n = 1, sd = 0.5) #sqrt(con$prop[j + 1]))
			}
		}
		 else{
			  # ANC STATES
			if (gibbs_sampler==F){
				k <- j - 1
				y.prime <- y
				y.prime[k] <- y[k] + rnorm(n = 1, sd = 0.5) #sqrt(con$prop[j + 1]))	
			} 
			else {
				if (i>100 && useDPP==T){
					for (i_temp in dim(D)[1]:1){
						dpp_list = DDP_gibbs_sampler4(tree,sig2.prime,ind_sig2,x, a.prime, y.prime,D,prior_tbl,alpha_par_Dir,i_temp,1)
						sig2.prime    = dpp_list[[1]]
						ind_sig2      = dpp_list[[2]]
						alpha_par_Dir = dpp_list[[3]]

						dpp_list = DDP_gibbs_sampler4(tree,sig2.prime,ind_sig2,x, a.prime, y.prime,D,prior_tbl,alpha_par_Dir,i_temp,2)
						sig2.prime    = dpp_list[[1]]
						ind_sig2      = dpp_list[[2]]
						alpha_par_Dir = dpp_list[[3]]
					}
					
				}
				vector_tip_root_nodes_values = c(x, a.prime, y.prime)
				y_temp = runGibbs(sig2.prime[ind_sig2], vector_tip_root_nodes_values,D,prior_tbl,get_expected=0)
				a.prime= y_temp[1]
				y.prime= y_temp[-1]
				gibbs=1
			}

		}

	       # calc post
		L.prime <- newlnLike(tree, c(x, a.prime, y.prime), sig2.prime[ind_sig2],D)
	
		Pr.prime <- log.prior(sig2.prime[ind_sig2], a.prime, y.prime)
	
		if ( (sum(Pr.prime) + sum(L.prime) - sum(Pr) - sum(L) + hastings) >= log(runif(1,0,1)) || gibbs==1){    
			y =    y.prime
			L =    L.prime
			Pr =   Pr.prime
			a =    a.prime
			sig2 = sig2.prime
		}
 
	     if (i%%sample == 0) {
			rates_temp=sig2[ind_sig2]
			#   rel_err =  (c(a, y) - true_anc)
			cat(c(i,sum(L)+sum(Pr), sum(L),sum(Pr), mean(sig2[ind_sig2]), length(sig2), a, alpha_par_Dir, rates_temp[IND_edge], y, "\n"),sep="\t", file=logfile, append=T) 
	    }
    
}
}

################################## END   MCMC ###################################


############################### SIMULATE DATA #######################################
sim_data <- function(ntips=20,s2=0.1,variable_rate=0){
	print("Simulating trees")
	tree<-ladderize(pbtree(n=ntips, scale=1))
	print("Simulating data")
	mod_tree <- tree
	
	x=sample.int(2,size= length(tree$tip.label),replace=T)
	names(x)= tree$tip.label
	
	# Martha\s code
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
	col = sample(colours(),length(phy$edge.length)) 
	names(col)=c(1:length(phy$edge.length))
	#plotSimmap(ttree, colors=col)
	
	# preparing rates for the simualtion 
	sigmas2= rep(s2,length(phy$edge.length))
	#alter = 1:10
	#sigmas2[alter] = sigmas2[alter]*5
	alter = sample(1:198,80) #5:80 # 
	if (variable_rate==1){
		sigmas2[alter] = sigmas2[alter]*15		
	}else if (variable_rate==2){
		alter = 5:80
		sigmas2[alter] = sigmas2[alter]*15		
	}
	#alter = sample(100:length(phy$edge.length),20) # 12:50
	#sigmas2[alter] = sigmas2[alter]/15
	# end Martha\s code
	
	#TE=as.matrix(tree$edge[order(tree$edge[,2]),])
	#names_l=c()
	#for (i in 1:length(sigmas2)){
	#	names_l[i] =paste(TE[i,1], ",",TE[i,2],sep="")
	#}
	names(sigmas2)  = NULL # names(col)#  #names_l
	
	#plot(mod_tree)
	#edgelabels()
	#nodelabels()
	full_data <- sim.rates(ttree, sig2=sigmas2, anc=0, nsim=1, internal=T, plot=T)
	#edgelabels()
	#nodelabels()
	plot.phylo(tree, edge.width=sigmas2*3,show.tip.label = F)
	print(sigmas2)
	
	#alter_ind = tree$edge[alter,2]
	
	
	#phenogram(ttree, full_data,add=F, col=col)
	#full_data<-fastBM(tree, sig2=s2, a=0, internal=T);
	
	#tree$edge = tree$edge*100
	
	return(c(tree,full_data,alter))
}

get_time <- function(){
	return(sum(as.numeric(strsplit(format(Sys.time(), "%X"),":")[[1]])*c(3600,60,1)))
}

start_MCMC_sim <- function(w_dir,sim_n,sig2,ntips,ngenerations,sampling_f,nGibbs_gen,Gibbs_sample,root_calibration = c(0,100), plot_res = F){
	setwd(w_dir)
	S = sim_data(ntips=ntips,s2=sig2,variable_rate=var_rate)
	tree= S[[1]]
	full_data= S[[2]]
	data <- S[[2]][names(S[[2]]) %in% tree$tip.label]
	true_anc = S[[2]][names(S[[2]]) %in% tree$edge]
	D <- build_table(tree, ntips,data)
	BR= branching.times(tree) # sorted from root to most recent
	dist_from_root = max(BR)-BR
	#phenogram(tree, full_data,add=F, col=rep("black",length(full_data)))
	alter_ind= S[[3]]	

	### LOG rel err
	prior_tbl = get_calibration_tbl(D,root_calibration)
	
	print("Starting MCMC...")

	t1 = get_time()
	logfile3 = sprintf("sim_%s_s2_%s_n_%s.log", sim_n, sig2, ntips)
	mcmc.gibbs4(tree, data,alter_ind, D, prior_tbl, true_anc, ngen= nGibbs_gen,gibbs_sampler=T,useVCV=F,logfile=logfile3,update_sig_freq=0.5,sample=Gibbs_sample,useDPP=T,dynamicPlot= T)
	cat("\nTime Gibbs:", get_time() - t1, sep="\t")

	if (plot_res==T){
		plot_traitgram <-function(logfile,add=F,col="black"){
			out_tbl = read.table(logfile,header=T)
			post_anc_states=c()
			for (i in 6:dim(out_tbl)[2]){
				 val = mean(out_tbl[250:dim(out_tbl)[1],i])
				 names(val) = length(post_anc_states)+length(data)+1
				 post_anc_states=append(post_anc_states,val)
			}
			print (post_anc_states[1])
			phenogram(tree, c(data, post_anc_states),ylim=c(-2,2),add=add,col=col)	
		}

		#phenogram(tree, full_data,ylim=c(-2,2),add=F)	
		plot_traitgram(logfile1,add=F,col="grey50")
		plot_traitgram(logfile2,add=T,col="red")
		plot_traitgram(logfile3,add=T,col="blue")
	}
}

start_MCMC_data <- function(w_dir,tree_file,data_file,prior_file,ngenerations,sampling_f,nGibbs_gen,Gibbs_sample,root_calibration = c(0,100), plot_res = F){
	setwd(w_dir)
	tree= read.tree(tree_file)
	ntips = length(tree$tip.label)
	raw_data= read.table(data_file,head=T)
	data = raw_data[,2]
	names(data) = raw_data[,1]
	D <- build_table(tree, ntips,data)
	BR= branching.times(tree) # sorted from root to most recent
	dist_from_root = max(BR)-BR

	### LOG rel err
	prior_tbl = read.table(prior_file) #get_calibration_tbl(D,root_calibration)
	prior_tbl = get_calibration_tbl(D,root_calibration)	
	
	print("Starting MCMC...")
	t1 = get_time()
	
	logfile3 = sprintf("%s_mcmc.log", strsplit(tree_file,".",fixed =T)[1])
	mcmc.gibbs4(tree, as.vector(data),alter_ind, D, prior_tbl, ngen= nGibbs_gen,gibbs_sampler=T,useVCV=F,logfile=logfile3,update_sig_freq=0.5,sample=Gibbs_sample,useDPP=T,dynamicPlot= T)
	cat("\nRun time:", get_time() - t1, sep="\t")

}




## GLOBAL
w_dir         = as.character(arg[1])   # working dir
sim_n         = as.integer(arg[2])     # simulation number
sig2          = as.double(arg[3])      # BM rate
ntips         = as.integer(arg[4])     # number of tips
ngen          = as.integer(arg[5])     # number of mcmc gen
sample        = as.integer(arg[6])     # sampling freq
nGibbs_gen    = as.integer(arg[7])     # number of mcmc gen Gibbs
Gibbs_sample  = as.integer(arg[8])     # sampling freq Gibbs

#plot_res= F
#root_calibration =c(0,100)

start_MCMC_sim(w_dir, sim_n, sig2, ntips,ngen,sample,nGibbs_gen,Gibbs_sample)

#___	start_MCMC("/Users/daniele/Documents/projects/fossilBM/", 1, 0.2, ntips=100,ngenerations=1000000,sampling_f=2000,nGibbs_gen=250000,Gibbs_sample=500)
#___	"/Users/daniele/Documents/projects/fossilBM/" 1 0.2 100 1000000 2000 250000 100
#___	
#___	
#___	
#___	
w_dir            = "/Users/daniele/Documents/projects/fossilBM/"
sim_n            = 1
sig2             = 0.1
ntips            = 100
ngen             = 5000
nGibbs_gen       = 250000
Gibbs_sample     = 100
root_calibration = c(0,100)
var_rate         = 1
print_f          = 100
useDPP           = T
dynamicPlot      = T


# empirical
w_dir = "/Users/daniele/Dropbox-personal/Dropbox/fossilbm/data/"
data_file = "latitude_mean.txt"
tree_file  = "tree_used.tre"










