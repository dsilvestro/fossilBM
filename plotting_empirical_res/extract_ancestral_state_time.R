
library(geiger)
library(phangorn)
library(phytools)
library(scales)

setwd("your_working_directory")
rda_files <- list.files(pattern=".rda")
log_files <- list.files(pattern=".log")


#### Function for generating the phenotypic values at branching times, produces and array of matrices with 
#### time, min, max, sd of phenotype
 
generating_array <- function (slices, phy, data_phen, anc_mcmc_sample) { 

  ##########################################################    
  ## Loop over mcmc steps 
  ##########################################################
  
  #results_all <- array(dim=c(length(slices), 4, reps))

  for (r in c(1:reps)){
    
  print(paste("mcmc step = ", r, sep=""))
	

	# trait sample mcmc r 
  x <- unlist ( c(data_phen, anc_mcmc_sample[r,]) )
    

  ###############################################################  
  ### Loop inside a mcmc step, over branching times, here slices 
  ###############################################################
  
      for(s in 1:length(slices)){
    
		  #print(paste("mcmc=", r, "slice=", s, sep=" "))
			t=slices[s]
      
		    #################################################
		  	### for the root slice, because trait=root state
		    #################################################
		  
			if(s == 1) {
		    
		  #  print(x[root])
			cat(x[root],t,"\n",sep="\t",file="data_trait_time.log",append=T)
		   
			} else { 
			  
		    #################################################################
		    ## for all the other steps, where slice summarized all branches
		    ################################################################
		  		
		    node_at_t <- phy$edge[ which(H[,1] == t)[1] , 1] 		
		    edges_at_t <-which(H[,2]< t & H[,1] > t)
		
		    ## taking all the phenotypes associated with node and edges 
		    ## the one at THE node slice
				
		    p_at_node <- as.numeric(x[node_at_t] ) 
		    
		    ## The others by bm calculation from node
		
		    b=phy$edge[edges_at_t,]
			
			# nodes that belongs to the branches cut at point t
			# trait values at nodes connecting branches x and y (line 1 in b)
			# collecting traits at this point time, but weigthed by bl from/to each node
        
		    ##############################################################################
		    # If only the node is crossing that time, for some trees is the case (length b=0)  
		    ##############################################################################
		    
			if (length(b)==0) {
		    
		   #  print(paste(c(p_at_node, "time", t)))  
			cat(p_at_node,t,"\n",sep="\t",file="data_trait_time.log",append=T)
		    					}
		    
		    ##############################################################################
		    # If only a branch crosses that time (length b=2), don't loop over rows of b  
		    ##############################################################################
		    
			if (length(b)==2) {
                
				trait= x [b]  # trait at the starting node and ending 
				# time passed from node until cutting point    
					
				p1 = (as.numeric(slices[which(names(slices)==(names(trait)[1]))])) - t 
				
				# need the time of second node or fossil p2 
				 
				# if fossil, need to take the edge heigth
				if (names(trait)[2] %in% phy$tip.label ) {
					
						p2 = t - H [ which (phy$edge[,2] == (which( phy$tip.label == names(trait)[2]) ) ) , 2]
							}	else {		
						p2 = t - (as.numeric(slices[which(names(slices)==(names(trait)[2]))])) 
							}
			
				wt = trait[1]+(trait[2]-trait[1]) * p1 / (p1+p2)	
				
				#t1 = p1 * trait[1]
				#t2 = p2 * trait[2]
				#wt = (t1+t2)/(p1+p2)
			
				#print(paste(c(wt, "time", t)))  
				cat(wt,t," \n",sep="\t",file="data_trait_time.log",append=T)
								
					} 
		    
			##############################################################################
			# if many branches crosses that time, loop over rows of b 
 			#############################################################################
			
			if (length(b) > 3) { 
				
					
				for (l in 1:nrow(b)) {
			
							trait= x[b[l,]]
							# time passed from node until cutting point    
							p1 = (as.numeric(slices[which(names(slices)==(names(trait)[1]))])) - t 
							# need the time of second node or fossil p2 
							# if fossil, need to take the edge heigth
							if (names(trait)[2] %in% phy$tip.label ) {
								p2 = t - H [ which (phy$edge[,2] == (which( phy$tip.label == names(trait)[2]) ) ) , 2]
									}	else {		
								p2 = t - (as.numeric(slices[which(names(slices)==(names(trait)[2]))])) 
									}
							
							wt = trait[1]+(trait[2]-trait[1]) * p1 / (p1+p2)	
							
							#t1 = p1 * trait[1]
							#t2 = p2 * trait[2]
							#wt = (t1+t2)/(p1+p2)
				 
						#	print(paste(wt,"many branches",t))
						
						    cat(wt,t,"\n",sep="\t",file="data_trait_time.log",append=T)	
							
							}
  
					}
		    
		    }
		
  
		    #############################
      } # close the slices loop  
		    #############################
		    
	
    #################################
    }   # close the mcmc steps loop
    #################################
  
	print("close mcmc")
  
}  # close function 


#############################################
###### Loop over files ######

cat("trait","time\n",sep="\t",file="data_trait_time.log", append=T) 


for (f in 1:length(rda_files)) {

	  # load original data

	  # robj_file = sprintf("%s/sim_%s.rda", wd, sim)
	  # robj_file = "platyrrhines0221corrected_names.trees1LAT.rda"
	  robj_file = rda_files[f]
	  print(robj_file)
	  load(robj_file)

	  # load MCMC samples
	  # logfile3= sprintf("%s/sim_%s.log",wd,sim)
	  # logfile3= "platyrrhines0221corrected_names.trees1LAT.log"
	  logfile3= log_files[f]
	  # Summarizing sigma2 values 			
	  out_tbl = read.table(logfile3,header=T)
	  burnin= round(0.25*dim(out_tbl)[1])
	  ind_anc_col = grep('anc_', colnames(out_tbl), value=F)
	  anc_mcmc = out_tbl[burnin:dim(out_tbl)[1], ind_anc_col]

	  phy=tree
	  # node.heights
	  root<-length(phy$tip)+1
	  H<-nodeHeights(phy)
	  H <- -(H-max(H))
  
	  # slices of time to loop over 
	  slices <- as.numeric(H[,1])
	  names(slices)<- paste('anc_', phy$edge[,1], sep="")
	  slices <- sort(slices[unique(names(slices))], decreasing = T)

	  ## the number of samples from the mcmc , or we can take the 95% of the mcmc after burnin 
	  reps=20
	  anc_mcmc_sample = anc_mcmc [ sample(1:dim(anc_mcmc)[1], size=reps, replace=F) , ] 
	  data_phen <- data[phy$tip.label]

	  ####### Running the function over a tree with 100 mcmc steps ######################
	  ####### This produces an array with matrices that have time, min, max, 
	
      results_all_test <- generating_array (slices, phy, data_phen, anc_mcmc_sample)

		
	  ###################################################################################
  
	 # rangerepf <- rangeForpolygon (results_all_test)  
  
	#  assign( paste("range", f, sep="_"), rangerepf)
}

