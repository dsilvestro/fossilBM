require(phytools);
library(diversitree)

update_multiplier_proposal <- function(i,d){
	u = runif(1,0,1)
	l = 2*log(d)
	m = exp(l*(u-0.5))
 	ii = i * m
	U=log(m)
	return(c(ii, U))
}


fN<-function(xa, xb, vpa, vpb, sigma2) {
  #the normal density. Sigma2 is the variance. Same as log(dnorm(xa-xb, 0, sqrt(sigma * sigma * (vpa + vpb))))
  return( dnorm((xa-xb),mean=0,sd= sqrt(sigma2 * (vpa + vpb)), log=T) )
  #return((-(xa - xb)^2 / (2 * sigma2 * (vpa + vpb)) - log(sqrt(2 * pi * sigma2 * (vpa + vpb)))));
}

fN2<-function(xa, xb, vpa, vpb, sigma2, anc) {
  #same as fN but with a known ancestor instead of xa-xb, see Joe's book eqn 23.10 (some mistakes in the denominator though in his equation)
  return(dnorm(xa, mean=anc, sd= sqrt(vpa*sigma2),log=T) + dnorm(xb, mean=anc, sd= sqrt(vpb*sigma2),log=T)  )
  #return(-(((xa - anc)^2 / vpa) + ((xb - anc)^2 / vpb)) / (2 * sigma2) - log(sqrt(2 * pi * sigma2 * (vpa + vpb))));
}


newlnLike<-function(tree, vector_tip_root_nodes_values, sigma2) {
  vec_values = vector_tip_root_nodes_values

  #conditional likelihood vectors
  condvec<-rep(0,length(vec_values))

  # loop over Johnatan's tbl (REVERSE ORDER)
  for (i in dim(D)[1]:1){
	anc_ind <- D[i,1];
	a_ind   <- D[i,2];  # index of descendants
	b_ind   <- D[i,3];  # index of descendants
	vpa     <- D[i,4];  # br length
	vpb     <- D[i,5];  # br length
	anc = vec_values[anc_ind]
	a = vec_values[a_ind]
	b = vec_values[b_ind]
	condvec[anc_ind] = fN2(a,b, vpa+S_vec[a_ind], vpb+S_vec[b_ind], sigma2, anc) # + condvec[a_ind] + condvec[b_ind];
	#print(c(i, anc_ind))
	
   }
   lik = sum(condvec) + fN(anc,anc, S_vec[anc_ind], 0, sigma2) 
   #print(condvec)
   return(lik) 
  
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


get_S_vec <- function (D){
	S =rep(0,(dim(D)[1]+dim(D)[1]+1))
	
        for (i in dim(D)[1]:1){
	      	anc_ind <- D[i,1];
	      	a_ind   <- D[i,2];  # index of descendants
	      	b_ind   <- D[i,3];  # index of descendants
	      	vpa     <- D[i,4];  # br length
	      	vpb     <- D[i,5];  # br length
	      	
		
		S[anc_ind] = ((vpa + S[a_ind]) * (vpb + S[b_ind])) / ((vpa + S[a_ind]) + (vpb + S[b_ind]))
		
	    #  	print(S)
	}
	return(S)
}


calc_prior_dist_root_calibration <- function (D,sigma2){
	trend=0
	mu_anc = root_calibration[1] + dist_from_root *trend
	s2_anc = root_calibration[2] + dist_from_root *sigma2
	prior_tbl = cbind(mu_anc,s2_anc)
	return(prior_tbl)
}


################ RUN GIBBS FUNCTIONS
get_joint_mu_s2 <- function (mu_f,s2_f,mu_g,s2_g){
	s2_fg = (s2_f*s2_g)/(s2_f+s2_g)
	mu_fg = (mu_f/s2_f + mu_g/s2_g) * s2_fg
	return(c(mu_fg,s2_fg))
}

runGibbs <- function(sigma2, vector_tip_root_nodes_values) {
	
	prior_tbl = calc_prior_dist_root_calibration(D,sigma2)
	prior_tbl[2:dim(prior_tbl)[1],1]=0
	prior_tbl[2:dim(prior_tbl)[1],2]=100	
	
	vec_values = vector_tip_root_nodes_values
        # loop over Johnatan's tbl from most recent to root
        for (i in dim(D)[1]:1){
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

		desc_lik = get_joint_mu_s2(a, (vpa+S_vec[a_ind])*sigma2, b, (vpb+S_vec[b_ind])*sigma2)
		
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
		
			stem_val = vec_values[stem_ind]		
			stem_lik = c(stem_val, (stem_brl+S_vec[stem_ind])*sigma2)
			#stem_lik = get_joint_mu_s2(prior_tbl[i-1,][1],prior_tbl[i-1,][2], stem_lik[1],stem_lik[2] )
			#get_joint_mu_s2(calibrated_prior_mu,calibrated_prior_s2,stem_val,(stem_brl+S_vec[stem_ind]) *sigma2)
			
			lik_val = get_joint_mu_s2(stem_lik[1],stem_lik[2],desc_lik[1],desc_lik[2])
			prior_prm = prior_tbl[which(rownames(prior_tbl)==anc_ind),]
		}
		
		else { # ROOT
			lik_val = desc_lik
			mu_0= anc
			s2_0= S_vec[anc_ind]
			prior_prm = prior_tbl[which(rownames(prior_tbl)==anc_ind),]	
			prior_prm = get_joint_mu_s2(mu_0,s2_0,prior_prm[1],prior_prm[2])			
		}
		
		
		
		#if (i==1){
		#	prior_prm = root_calibration
		#}else{
		#	prior_prm = c(0,100)
		#}
		
		
		#print(c(i,prior_prm[1],prior_prm[2]))
		
		if (i>=1){ 
			post_prm= get_joint_mu_s2(lik_val[1],lik_val[2], prior_prm[1], prior_prm[2])
			vec_values[anc_ind] = rnorm(1, mean=post_prm[1], sd=sqrt(post_prm[2]))
		}
		# else{
		# 	vec_values[anc_ind] = rnorm(1, mean=prior_prm[1], sd=sqrt(prior_prm[2]))
		# }
         }
         #print(condvec)
         return(vec_values[D[,1]] )
}



################################## START MCMC ###################################
######################################################################
mcmc.gibbs4 <- function (tree, x, ngen = 100000, control = list(), gibbs_sampler=T,useVCV=F, sample=100,logfile="log",update_sig_freq=0.5) 
{
    cat(c("it", "posterior","likelihood","prior", "sig2", "root", D[-1,1],"\n"),sep="\t", file=logfile, append=F)
    
    temp <- phyl.vcv(as.matrix(x), vcv(tree), 1)
    sig2 <- 0.5 #temp$R[1, 1]
    a <- temp$alpha
    y <- rep(a, tree$Nnode - 1)
    
    # lik function
    if (useVCV ==T){
	    likelihood <- function(C, invC, detC, x, sig2, a, y) {
	        z <- c(x, y) - a
	        logLik <- -z %*% invC %*% z/(2 * sig2) - nrow(C) * log(2 * 
	            pi)/2 - nrow(C) * log(sig2)/2 - detC/2
	        return(logLik)
	    }    	
    }
    else{
	    likelihood <- function(C, invC, detC, x, sig2, a, y) {
		    vector_tip_root_nodes_values = c(x, a, y)
		newlnLike(tree, vector_tip_root_nodes_values, sig2)
	    }    	
    }

    # priors
    log.prior <- function(sig2, a, y) {
	prior_sig2 <- dexp(sig2, rate = 1/100, log = TRUE) 
	# + sum(dnorm(c(a, y), mean = pr.mean[2:length(pr.mean)], sd = sqrt(pr.var[1 + 1:tree$Nnode]), log = TRUE))
	prior_tbl = calc_prior_dist_root_calibration(D,sig2)
	prior_root = sum(dnorm(c(a), mean = prior_tbl[,1], sd = prior_tbl[,2], log = T))
	prior_anc = sum(dnorm(c(a, y), mean = 0, sd = 100, log = T))
    	return(prior_sig2+prior_root+prior_anc)
    }
	
    C <- vcvPhylo(tree)
    if (any(tree$edge.length <= (10 * .Machine$double.eps))) 
        stop("some branch lengths are 0 or nearly zero")
    invC <- solve(C)
    detC <- determinant(C, logarithm = TRUE)$modulus[1]
    x <- x[tree$tip.label]
    if (is.null(names(y))) 
        names(y) <- length(tree$tip) + 2:tree$Nnode
    else y[as.character(length(tree$tip) + 2:tree$Nnode)]
    L <- likelihood(C, invC, detC, x, sig2, a, y)
    Pr <- log.prior(sig2, a, y)
    
    # START MCMC
    for (i in 1:ngen) {
	    
        y.prime     = y
        L.prime     = L
        Pr.prime    = Pr
        a.prime     = a
        sig2.prime  = sig2
        gibbs=0
	hastings=0
	    
	if (i%%1000 == 0) {print (round(c(i,L,sig2),2))}
	    
        j <- (i - 1)%%(tree$Nnode + 1)
        #if (j == 0) {
	rr= runif(1,0,1)
	if (rr<update_sig_freq) {
		#sig2.prime <-  abs(sig2 + rnorm(n = 1, sd = sqrt(con$prop[j + 1]))) 
		sig2_update <-  update_multiplier_proposal(sig2,1.2)
		sig2.prime = sig2_update[1]
		hastings = sig2_update[2]
    	}
         else {   # ANC STATES
		if (gibbs_sampler==F){
			a.prime <- a + rnorm(n = 1, sd = 0.5) #sqrt(con$prop[j + 1]))
			k <- j - 1
			y.prime <- y
			y.prime[k] <- y[k] + rnorm(n = 1, sd = 0.5) #sqrt(con$prop[j + 1]))	
		} 
		else {
			vector_tip_root_nodes_values = c(x, a, y)
			y_temp = runGibbs(sig2, vector_tip_root_nodes_values)
			a.prime= y_temp[1]
			y.prime= y_temp[-1]
			gibbs=1			
		}

       }

       # calc post

	
	L.prime <- likelihood(C, invC, detC, x, sig2.prime, a.prime,  y.prime)
	Pr.prime <- log.prior(sig2.prime, a.prime, y.prime)
	
	if ( (Pr.prime + L.prime - Pr - L + hastings) >= log(runif(1,0,1)) || gibbs==1){    
	#if (2>1){
		y =    y.prime
		L =    L.prime
		Pr =   Pr.prime
		a =    a.prime
		sig2 = sig2.prime
	}
 
 
     if (i%%sample == 0) {
	    cat(c(i,L+Pr, L,Pr, sig2, a, y,"\n"),sep="\t", file=logfile, append=T)
 
    }
    
}
message("Done.")
}



################################## END   MCMC ###################################

### LOG rel err
# sig2 in file name
# n tips in file name
# method in file name
# add ags for sig2, tips



tree<-pbtree(n=20, scale=1);
full_data<-fastBM(tree, sig2=0.5, a=0, internal=T);
data <- full_data[names(full_data) %in% tree$tip.label]
true_anc = full_data[names(full_data) %in% tree$edge]
D <- build_table(tree, length(tree$tip.label),data)
S_vec <- get_S_vec(D) # vector with extra-variances
BR= branching.times(tree) # sorted from root to most recent
dist_from_root = max(BR)-BR
root_calibration =c(0,100)


# run VCV likelihood
logfile1 = "/Users/daniele/Desktop/logVCV.log"
mcmc.gibbs4(tree, data, 25000,gibbs_sampler=F,useVCV=T,logfile=logfile1,update_sig_freq=0.1)

# run MH likelihood
logfile2 = "/Users/daniele/Desktop/logMH.log"
mcmc.gibbs4(tree, data, 25000,gibbs_sampler=F,useVCV=F,logfile=logfile2,update_sig_freq=0.1)

# run Gibbs likelihood
logfile3 = "/Users/daniele/Desktop/logGibbs.log"
mcmc.gibbs4(tree, data, 25000,gibbs_sampler=T,useVCV=F,logfile=logfile3,update_sig_freq=0.75)

plot_traitgram <-function(logfile,add=F,col="black"){
	out_tbl = read.table(logfile,header=T)
	post_anc_states=c()
	for (i in 6:dim(out_tbl)[2]){
		 val = mean(out_tbl[50:dim(out_tbl)[1],i])
		 names(val) = length(post_anc_states)+length(data)+1
		 post_anc_states=append(post_anc_states,val)
	}
	print (post_anc_states[1])
	phenogram(tree, c(data, post_anc_states),ylim=c(-2,2),add=add,col=col)	
}

phenogram(tree, full_data,ylim=c(-2,2),add=F)	

plot_traitgram(logfile1,add=T,col="grey50")
plot_traitgram(logfile2,add=T,col="red")
plot_traitgram(logfile3,add=T,col="blue")

