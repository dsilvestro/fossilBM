#!/usr/bin/Rscript
# version 20150713


# module add R/3.0.2
# cd /scratch/ul/monthly/dsilvest/BM && RScript mcmcFossilBM-BD.R --wd /scratch/ul/monthly/dsilvest/BM --j %s --v 8 --s 2 --r 0.1 --t 100


arg <- commandArgs(trailingOnly=TRUE)

require(phytools)
require(diversitree)
#library(gtools)
library(optparse)
##################################################

option_list <- list(                   
	make_option("--n",  type="integer",   default=250000,help=("MCMC iterations (%default)"),metavar="iterations"),
	make_option("--t",  type="integer",   default=50,help=("Tree size (%default)"),metavar="ntips"),
	make_option("--v",  type="double",    default=16,help=("Magnitude rate shift (%default)"),metavar="xfold"),    
	make_option("--r",  type="double",    default=0.1,help=("Baseline rate (%default)"),metavar="rate"),    
	make_option("--j",  type="integer",   default=1,help=("Replicate (%default)"),metavar="replicate"),
	make_option("--wd", type="character", default="",help=("working dir"),metavar="wd"),
	make_option("--s",  type="integer",   default=1,help=("no. shifts (%default)"),metavar="shifts")
)

parser_object <- OptionParser(usage = "Usage: %prog [Options]", option_list=option_list, description="...")
opt <- parse_args(parser_object, args = commandArgs(trailingOnly = TRUE), positional_arguments=TRUE)
	
ntips            = opt$options$t
nGibbs_gen       = opt$options$n
Gibbs_sample     = 500
root_calibration = c(0,100)
print_f          = 250
xfold            = opt$options$v
dynamicPlot      = F
w_dir            = opt$options$wd #"/Users/daniele/Dropbox-personal/Dropbox/fossilizedBM"
rep              = opt$options$j
sim_n            = 1
sig2             = opt$options$r
n_shifts         = opt$options$s

### TEMP SETTINGS
ntips            = 25
nGibbs_gen       = 10000
Gibbs_sample     = 50
root_calibration = c(0,0.05)
print_f          = 500
xfold            = 1
dynamicPlot      = T
w_dir            = "/Users/daniele/Dropbox-personal/Dropbox/fossilizedBM"
rep              = 1
sim_n            = 1
sig2             = 0.25
mu0              = 0.5
n_shifts         = 0




setwd(w_dir)

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

newlnLike<-function(tree, vector_tip_root_nodes_values, sigma2,D,mu0) {
	#S_vec=rep(0,10000)
        #vec_values = vector_tip_root_nodes_values
        #
        ##conditional likelihood vectors
        #condvec<-rep(0,length(vec_values))
        #
        ## loop over Johnatan's tbl (REVERSE ORDER)
        #for (i in dim(D)[1]:1){
        #	anc_ind <- D[i,1];
        #	a_ind   <- D[i,2];  # index of descendants
        #	b_ind   <- D[i,3];  # index of descendants
        #	vpa     <- D[i,4];  # br length
        #	vpb     <- D[i,5];  # br length
        #	anc = vec_values[anc_ind]
        #	a = vec_values[a_ind]
        #	b = vec_values[b_ind]
	#	print (c(a,b, vpa+S_vec[a_ind], vpb+S_vec[b_ind], sigma2, anc))
        #	condvec[anc_ind] = fN2(a,b, vpa+S_vec[a_ind], vpb+S_vec[b_ind], sigma2, anc) # + condvec[a_ind] + condvec[b_ind];
        #	print(c(i, anc_ind,condvec[anc_ind]))
        # }
        # lik = condvec #+ fN(anc,anc, S_vec[anc_ind], 0, sigma2) 
        #
	
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

runGibbs <- function(sigma2, vector_tip_root_nodes_values,D,prior_tbl,mu0,get_expected=0) {
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

		
		#desc_lik = get_joint_mu_s2(a+vpa, vpa*sigma2[sig_a_ind], b, (vpb)*sigma2[sig_b_ind])

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
			
			#if (stem_ind>ntips){
			sig_stem_ind=stem_ind-1
			#}else{sig_stem_ind=stem_ind}
			
			
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
	#prior_sig2 <- sum(dcauchy(sig2, loc=0,scale=1, log = TRUE) )
	prior_sig2 = sum(dexp(sig2, 0.5, log = TRUE) ) #sum(dgamma(sig2, 2,1, log = TRUE) )
	prior_root = sum(dnorm(c(a), mean = prior_tbl[,1], sd = prior_tbl[,2], log = T))
	prior_anc  = sum(dnorm(c(a, y), mean = 0, sd = 100, log = T))
	prior_mu0  = sum(dnorm(mu0, mean = 0, sd = 10, log = T))
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

calc_deathRates <- function (tree, trait_states, sig2, ind_sig, mu0, ind_mu0,D,IND_edge,new_ind,desc_index_list){
	list_deathrates=c(0) # can't kill first
	current_lik = sum(newlnLike(tree, trait_states, sig2[ind_sig],D,mu0[ind_mu0]))
	#print(length(current_lik))
	#__ par(mfrow=c(1,length(unique(ind_sig))-1))
	for (i in 2:length(unique(ind_sig))){
		multi_ind = ind_sig
		multi_ind_edge  = multi_ind[IND_edge]
		# get edges with index i
		Yi = which(multi_ind_edge==i)
		#__ cat("\ndead edge:",Y,"ind:" ,i,"\n")
		crown_node_ind_i = min(tree$edge[Yi,1])
		#__ cat("\nind crown:",crown_node_ind_i,"\n")
		if (crown_node_ind_i==(ntips+1)){ # ROOT
			ind_sig_of_ancestor=1
		} else{
			ind_sig_of_ancestor= multi_ind_edge[which(tree$edge[,2]==crown_node_ind_i)]
		}

		new_ind_sig = ind_sig
		new_ind_sig[new_ind_sig==i] = ind_sig_of_ancestor
		
		rates_temp= sig2[new_ind_sig]
		#__ print(c(i, sig2[i], sig2[ind_sig_of_ancestor]) )
		#__ plot.phylo(tree, edge.width=rates_temp[IND_edge],show.tip.label = F)
		
		#y_temp = runGibbs(sig2[ind_sig], trait_states,D,prior_tbl,get_expected=0)
		#a.prime= y_temp[1]
		#y.prime= y_temp[-1]
		#trait_states_temp = trait_states
		#trait_states_temp[(ntips+1):length(trait_states)] = c(a.prime,y.prime)			
		
		list_deathrates[i] = sum(newlnLike(tree, trait_states, rates_temp,D,mu0[ind_mu0]))
	}
	list_deathrates = exp(list_deathrates - current_lik)
	list_deathrates[1] = 0
	return(list_deathrates)
}

run_BDMCMC <- function(tree, trait_states, sig2, ind_sig, mu0, ind_mu0,D,IND_edge,new_ind,desc_index_list){
	cont_time=0
	len_cont_time=1
	birthRate=1
	while (cont_time<len_cont_time){
		#if (cont_time>=0){ # run gibbs update of anc states after each BD event
		#	y_temp = runGibbs(sig2[ind_sig], trait_states,D,prior_tbl,get_expected=0)
		#	a.prime= y_temp[1]
		#	y.prime= y_temp[-1]
		#	trait_states[(ntips+1):length(trait_states)] = c(a.prime,y.prime)			
		#}
		
		
		deathRate = rep(0, length(unique(ind_sig)))
		n_likBD   = rep(0, length(unique(ind_sig)))
		
		#print( length(unique(ind_sig)))
		
		if (length(unique(ind_sig))>1){
			deathRate = calc_deathRates(tree, trait_states, sig2, ind_sig, mu0, ind_mu0, D,IND_edge,new_ind,desc_index_list)
		}
		
		deltaRate=sum(deathRate)
		
		cont_time = cont_time + rexp(1, rate = min((deltaRate+birthRate), 100000))
		if (cont_time>len_cont_time){
			"pass"
		}else{
			Pr_birth= birthRate/(birthRate+deltaRate)
			Pr_death= 1-Pr_birth
			
			#__ cat("\nB:",Pr_birth, "D:",Pr_death)
			#__ cat("\nDELTA", deathRate, "\n")
			
			if (runif(1)<Pr_birth){ # ADD PARAMETER
				rates_temp= sig2[ind_sig]
				#__ plot.phylo(tree, edge.width=rates_temp[IND_edge],show.tip.label = F)
				n_ind  = sample(1:length(ind_sig),1)
				
				# sample an index that will be split into two indexes
				#tab = table(ind_sig)
				#n_ind=sample(as.numeric(names(tab[tab>1])),1) # all indexes have = chance of being split if the contain > 1 branch
				#__ cat("\nchosen index", n_ind) #, "table", tab)
				
				
				multi_ind = ind_sig
				multi_ind_edge  = multi_ind[IND_edge]
				multi_ind_new = multi_ind_edge
				# index that is being affected
				split_ind = multi_ind_new[n_ind]
				
				ind_desc = desc_index_list[[n_ind]]
				for (j in ind_desc){
					#print(c(j,multi_ind_new[j]))
					if (multi_ind_new[j]==split_ind){multi_ind_new[j] = max(ind_sig)+1}
				}
				cat()
				#___ # alternative - NO LOOP
				#___ multi_ind_new[which(ind_desc==split_ind)] =  max(ind_sig)+1
				

				# BDMCMC PROPOSAL
				sig2=c(sig2,rexp(1,0.5)) # sample from rate prior
				
				#temp_r = update_multiplier_proposal(sig2[1],10) # rate as multiplier update of baseline rate
				#sig2=c(sig2,temp_r[1]) 
				
				
				#__ rates_temp1 = sig2[multi_ind_new]
				#__ plot.phylo(tree, edge.width=rates_temp1,show.tip.label = F)
				
				# ind_sig in the original order
				ind_sig=multi_ind_new[new_ind]
				
				#__ cat("\npre", unique(multi_ind), " post",unique(ind_sig), " sig2", length(sig2) )
				
				if (length(sig2)>length(unique(ind_sig))){
					# the new index overlaps completely with the previous
					# remove unused sig2
					a = unique(multi_ind)
					b = unique(ind_sig)
					old_ind_that_disappeared = min(as.numeric(names(table(c(a,b))[table(c(a,b))==1])))
					sig2 = sig2[-old_ind_that_disappeared]
					# rescale indexes
					ind_sig[ind_sig>(old_ind_that_disappeared)]=ind_sig[ind_sig>(old_ind_that_disappeared)]-1
				}
				#__ cat(" FINAL:", unique(ind_sig), " sig2", length(sig2))
				
				# rates back to original order
				#__ rates_temp= sig2[ind_sig]
				#__ plot.phylo(tree, edge.width=rates_temp[IND_edge],show.tip.label = F)
				
			} 
			else{ # REMOVE PARAMETER
				i=random_choice_P(deathRate)[2]
				
				multi_ind = ind_sig
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
					cat("\nAHAHAHAHAHAHAHAHAH",length(unique(ind_sig)),max(ind_sig),length(deathRate),"unique",unique(ind_sig))
				}
				
				if (crown_node_ind_i==(ntips+1)){ # ROOT
					ind_sig_of_ancestor_i=1
				} else{
					ind_sig_of_ancestor_i= multi_ind_edge[which(tree$edge[,2]==crown_node_ind_i)]
				}
				
				# update indexes
				new_ind_sig = ind_sig
				new_ind_sig[new_ind_sig==i] = ind_sig_of_ancestor_i
				ind_sig = new_ind_sig
				
				# remove unused sig2
				sig2 = sig2[-i]
				
				# rescale indexes
				#__cat("\nproto", unique(multi_ind), " pre",unique(ind_sig) )
				ind_sig[ind_sig>i]=ind_sig[ind_sig>i]-1
				#__cat(" post",unique(ind_sig),"\n")
				
				#__ rates_temp= sig2[ind_sig]
				#__ plot.phylo(tree, edge.width=rates_temp[IND_edge],show.tip.label = F)
				

			}
			
		}
		
	}
		
	return (list(sig2,ind_sig))
}

################################## START MCMC ###################################
mcmc.gibbs4 <- function (tree, x, D, prior_tbl,true_rate,ngen = 100000, control = list(),useVCV=F, sample=100,logfile="log",update_sig_freq=0.5,dynamicPlot = F,useFixPart=F,bdmcmc_freq=0.75,useTrend=F){
	x=data
	TE=as.matrix(tree$edge)
	cat(c("it", "posterior","likelihood","prior", "sig2", "mu0", "MAPE","sdAPE","K", paste("sig2", 1:length(TE[,1]),sep="_"), paste("anc",D[,1],sep="_"),"\n"),sep="\t", file=logfile, append=F)
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
			if (runif(1,0,1)<bdmcmc_freq){
				trait_states=c(x, a.prime, y.prime)
				bdmcmc_list   = run_BDMCMC(tree, trait_states, sig2.prime,ind_sig2, mu0.prime, ind_mu0,D,IND_edge,new_ind,desc_index_list)
				sig2.prime    = bdmcmc_list[[1]]
				ind_sig2      = bdmcmc_list[[2]]
				gibbs=1
			}
			vector_tip_root_nodes_values = c(x, a.prime, y.prime)
			y_temp = runGibbs(sig2.prime[ind_sig2], vector_tip_root_nodes_values,D,prior_tbl,mu0.prime[ind_mu0],get_expected=0)
			a.prime= y_temp[1]
			y.prime= y_temp[-1]
			gibbs=1
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
			#   rel_err =  (c(a, y) - true_anc)
			#MAE = mean(abs(c(a, y)-true_anc))
			MAPE = mean(abs((rates_temp[IND_edge]-true_rate))/true_rate)
			sdAPE = sd(abs((rates_temp[IND_edge]-true_rate))/true_rate)
			cat(c(i,sum(L)+sum(Pr), sum(L),sum(Pr), mean(sig2[ind_sig2]),mean(mu0[ind_mu0]), MAPE, sdAPE, length(sig2), rates_temp[IND_edge], a, y, "\n"),sep="\t", file=logfile, append=T) 
	    }
    
}
}

################################## END   MCMC ###################################

############################### SIMULATE DATA #######################################
sim_data_bmt <- function(ntips=20,s2=0.1,m0=0,root_value=0){
	tree<-ladderize(pbtree(n=ntips, scale=1))
	D = build_table(tree,ntips,data=NULL)
	D = as.data.frame(D)
	D = cbind(D,rep(0,ntips-1))
	D = cbind(D,rep(0,ntips-1))
	D = cbind(D,rep(0,ntips-1))
	colnames(D)[6]="val0"
	colnames(D)[7]="val1"
	colnames(D)[8]="val2"
	
	tip_data=rep(0,ntips)
	# constant sig2
	sigmas2= rep(s2,length(tree$edge.length))
	
	for (node_ind in 1:length(D$ancestor)){
		node =	D$ancestor[node_ind]
		if (node_ind==1){
			anc_val = root_value
		}else{anc_val = D$val0[node_ind]}
		
		br_lengths=c(D$branch.length1[node_ind],D$branch.length2[node_ind])
		desc_val = anc_val + rnorm(n=2,mean=br_lengths*m0,sd=sqrt(br_lengths*s2))
		desc1_id = D$descendant1[node_ind]
		desc2_id = D$descendant2[node_ind]
		
		D$val1[node_ind]=desc_val[1]
		D$val2[node_ind]=desc_val[2]
		
		if (desc1_id>ntips){
			D$val0[which(D$ancestor==desc1_id)] = desc_val[1]
		}else{tip_data[desc1_id]= desc_val[1]}
		
		if (desc2_id>ntips){
			D$val0[which(D$ancestor==desc2_id)] = desc_val[2]
		}else{tip_data[desc2_id]= desc_val[2]}
		
	}
	all_data = c(tip_data, D$val0)
	all_names = c(tree$tip.label,D$ancestor)
	names(all_data) = all_names
	phenogram(tree, all_data,add=F, col="black", main="Traigram")
	return(c(tree,all_data,sigmas2,m0))
}

sim_data <- function(ntips=20,s2=0.1,xfold=16,n_shifts=1){
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
	#col = sample(colours(),length(phy$edge.length)) 
	col = rep("black",length(phy$edge.length))
	names(col)=c(1:length(phy$edge.length))
	#plotSimmap(ttree, colors=col)
	
	# preparing rates for the simualtion 
	sigmas2= rep(s2,length(phy$edge.length))
	
	
	#### HERE!!!
	#__ plot(tree)
	#__ nodelabels()
	#__ edgelabels()
	rnd_sample_tax = sample((1:ntips)+ntips, ntips)
	chosen_nodes =c()
	r_multipliers = runif(n_shifts,8,16)
	j=1
	if (n_shifts>=1){
		for (s in 1:n_shifts){
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
			print(c(i, length(intersect(chosen_nodes,ind_desc_nodes))))
			chosen_nodes= append(chosen_nodes,ind_desc_nodes)
			sigmas2[ind_desc_edges] = s2*r_multipliers[s] # (xfold*s)
		}
	}else{ind_desc_edges=c(0)}
	
	names(sigmas2)  = NULL # names(col)#  #names_l
	
	#plot(mod_tree)
	#edgelabels()
	#nodelabels()
	full_data <- sim.rates(ttree, sig2=sigmas2, anc=0, nsim=1, internal=T, plot=F)
	
	plot.phylo(tree, edge.width=sigmas2*3,show.tip.label = F, main="True rates")
	#print(sigmas2)
	
	#alter_ind = tree$edge[alter,2]
	
	
	edge_w=abs(full_data[ttree$edge[,1]] - full_data[ttree$edge[,2]])/sqrt(tree$edge.length)
	plot.phylo(tree, edge.width=edge_w*3,show.tip.label = F,main="Empirical rates")
	#full_data<-fastBM(tree, sig2=s2, a=0, internal=T);
	
	#tree$edge = tree$edge*100
	phenogram(ttree, full_data,add=F, col=col, main="Traigram")
	
	return(c(tree,full_data,sigmas2))
}

get_time <- function(){
	return(sum(as.numeric(strsplit(format(Sys.time(), "%X"),":")[[1]])*c(3600,60,1)))
}

start_MCMC_sim <- function(w_dir,sim_n,sig2,ntips,ngenerations,sampling_f,nGibbs_gen,Gibbs_sample,root_calibration = c(0,100), plot_res = F){
	setwd(w_dir)
	ntips =50
	S = sim_data(ntips=ntips,s2=0.1,xfold=1,n_shifts=0)
	
	## SIM WITH TREND
	S = sim_data_bmt(ntips=ntips,s2=0.1,m0=0.75,root_value=0)
	
	tree= S[[1]]
	full_data= S[[2]]
	data <- S[[2]][names(S[[2]]) %in% tree$tip.label]
	true_anc = S[[2]][names(S[[2]]) %in% tree$edge]
	D <- build_table(tree, ntips,data)
	BR= branching.times(tree) # sorted from root to most recent
	dist_from_root = max(BR)-BR
	#phenogram(tree, full_data,add=F, col=rep("black",length(full_data)))
	true_sigmas = S[[3]]

	### LOG rel err
	prior_tbl = get_calibration_tbl(D,root_calibration)
	
	logfile3 = sprintf("sim_%s_s2_%s_n_%sT.log", sim_n, sig2, ntips)
	mcmc.gibbs4(tree, data, D, prior_tbl, true_rate=true_sigmas, ngen= nGibbs_gen,useVCV=F,
		logfile=logfile3,update_sig_freq=0.5,sample=Gibbs_sample,dynamicPlot= F,bdmcmc_freq=0.5,useTrend=T)


	# Test BM simulation
	#S = sim_data_bmt(ntips=ntips,s2=0.25,m0=0,root_value=0)
	#tree= S[[1]]
	#full_data= S[[2]]
	#data <- S[[2]][names(S[[2]]) %in% tree$tip.label]	
	##library(geiger)
	#fitContinuous(tree, data, model = c("BM"))$opt$sigsq
	    
		
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
	mcmc.gibbs4(tree, as.vector(data),alter_ind, D, prior_tbl, ngen= nGibbs_gen,useVCV=F,logfile=logfile3,update_sig_freq=0.5,sample=Gibbs_sample,dynamicPlot= T)
	cat("\nRun time:", get_time() - t1, sep="\t")

}




#################### RUN SIMULATION
resfile = sprintf("sim_%s_s2_%s_n_%s_x_%s.pdf", rep, sig2, ntips, xfold)
pdf(file=resfile,width=0.9*15, height=1.05*15)
par(mfrow=c(3,3))
S = sim_data(ntips=ntips,s2=sig2,xfold=xfold,n_shifts=n_shifts)
tree= S[[1]]
full_data= S[[2]]
data <- S[[2]][names(S[[2]]) %in% tree$tip.label]
true_anc = S[[2]][names(S[[2]]) %in% tree$edge]
D <- build_table(tree, ntips,data)
BR= branching.times(tree) # sorted from root to most recent
dist_from_root = max(BR)-BR
#phenogram(tree, full_data,add=F, col=rep("black",length(full_data)))
alter_ind= S[[3]]
true_sigmas = S[[3]]

### LOG rel err
prior_tbl = get_calibration_tbl(D,root_calibration)

print("Starting MCMC...")
t1 = get_time()
logfile3 = sprintf("sim_%s_s2_%s_n_%s_x_%s.log", rep, sig2, ntips, xfold)
mcmc.gibbs4(tree, data, D, prior_tbl, true_rate=true_sigmas, ngen= nGibbs_gen,bdmcmc_freq=0.5,logfile=logfile3,update_sig_freq=0.5,
	sample=Gibbs_sample,dynamicPlot=dynamicPlot,useTrend=T) 
cat("\nTime Gibbs:", get_time() - t1, "\n", sep="\t")

logfile2 = sprintf("sim_%s_s2_%s_n_%s_x_%s_const.log", rep, sig2, ntips, xfold)
mcmc.gibbs4(tree, data, D, prior_tbl, true_rate=true_sigmas, ngen= nGibbs_gen/10,bdmcmc_freq=0  ,logfile=logfile2,update_sig_freq=0.5,
	sample=Gibbs_sample/10,dynamicPlot=dynamicPlot,useTrend=F)

	
# plot res
out_tbl = read.table(logfile3,header=T)
ind_sig2_col = grep('sig2_', colnames(out_tbl), value=F)
burnin= round(0.25*dim(out_tbl)[1])
rates_temp_median = apply(out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=median,2)
rates_temp_mean = apply(out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=mean,2)
MAPE_median = mean(abs((rates_temp_median-true_sigmas))/true_sigmas)
sdAPE_median = sd(abs((rates_temp_median-true_sigmas)))
MAPE_mean = mean(abs((rates_temp_mean-true_sigmas))/true_sigmas)
sdAPE_mean = sd(abs((rates_temp_mean-true_sigmas)))
estK =out_tbl$K[burnin:dim(out_tbl)[1]]
estK = round(quantile(estK,probs = c(0.025,0.975)))
plot.phylo(tree, edge.width=rates_temp_median*3, main=paste("Estimated rates",sprintf("(K: %s-%s)",estK[1],estK[2])),show.tip.label = F)
plot(true_sigmas,  main=paste("MAPE_avg:",round(MAPE_mean,3),"MAPE_med:",round(MAPE_median,3)),pch=21,col="darkblue",bg="blue",ylim=c(0,max(c(true_sigmas,rates_temp_median,rates_temp_mean))))
points(rates_temp_median,pch=21,col="darkred",bg="red")
points(rates_temp_mean,pch=23,col="#FFA500",bg="#FFA500")

ind_anc_col = grep('anc_', colnames(out_tbl), value=F)
anc_mean = apply(out_tbl[burnin:dim(out_tbl)[1], ind_anc_col],FUN=mean,2)
mape = round(mean(abs(anc_mean-true_anc)/(max(true_anc)-min(true_anc))),3)
plot(anc_mean, true_anc, xlab="estimated anc states", ylab="true anc states",main=paste("Anc states",sprintf("(MAPE: %s)",mape)))
abline(coef=c(0,1),lty=2)

out_tbl = read.table(logfile2,header=T)
ind_sig2_col = grep('sig2_', colnames(out_tbl), value=F)
burnin= round(0.25*dim(out_tbl)[1])
rates_temp_median = apply(out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=median,2)
rates_temp_mean = apply(out_tbl[burnin:dim(out_tbl)[1], ind_sig2_col],FUN=mean,2)
MAPE_median = mean(abs((rates_temp_median-true_sigmas))/true_sigmas)
sdAPE_median = sd(abs((rates_temp_median-true_sigmas)))
MAPE_mean = mean(abs((rates_temp_mean-true_sigmas))/true_sigmas)
sdAPE_mean = sd(abs((rates_temp_mean-true_sigmas)))
plot.phylo(tree, edge.width=rates_temp_median*3, main="Estimated rates",show.tip.label = F)
plot(true_sigmas,  main=paste("MAPE_avg:",round(MAPE_mean,3),"MAPE_med:",round(MAPE_median,3)),pch=21,col="darkblue",bg="blue",ylim=c(0,max(c(true_sigmas,rates_temp_median,rates_temp_mean))))
points(rates_temp_median,pch=21,col="darkred",bg="red")
points(rates_temp_mean,pch=23,col="#FFA500",bg="#FFA500")

ind_anc_col = grep('anc_', colnames(out_tbl), value=F)
anc_mean = apply(out_tbl[burnin:dim(out_tbl)[1], ind_anc_col],FUN=mean,2)
mape = round(mean(abs(anc_mean-true_anc)/(max(true_anc)-min(true_anc))),3)
plot(anc_mean, true_anc, xlab="estimated anc states", ylab="true anc states",main=paste("Anc states",sprintf("(MAPE: %s)",mape)))
abline(coef=c(0,1),lty=2)

n<-dev.off()


  
