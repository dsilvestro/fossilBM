require(phytools)

## Based on the function simBM written by Liam J. Revell 2011, 2013
## I slightly modified bits and pieces to make sure that it can handle two
##   vectors. Works both for ultrametric trees and phylograms.
simBM<-function(tree,a=0,mu=0,sig2=1,bounds=c(-Inf,Inf),internal=FALSE,nsim=1){
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
	x<-matrix(data=rnorm(n=ndraws,mean=rep(mu*nedges,nsim),sd=rep(sqrt(sig2*nedges),nsim)),nedges,nsim)
	
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

## small example

tree<-pbtree(n=20)

ratesMu<-fastBM(tree, mu=2, sig2=0.5, internal=TRUE)[-(length(tree$tip.label)+1)]
ratesSig2<-fastBM(tree, mu=1, sig2=0.1, internal=TRUE)[-(length(tree$tip.label)+1)]

simBM(tree, a=0, mu=ratesMu, sig2=ratesSig2, internal=TRUE)

ntree<-tree
ntree$edge.length<-rexp(length(tree$edge.length), 0.5)
simBM(ntree, a=0, mu=ratesMu, sig2=ratesSig2, internal=TRUE)