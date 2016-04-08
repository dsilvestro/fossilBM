setwd("~/Dropbox-personal/Dropbox/fossilizedBM/FBD_simulator")
library(phytools)

tree<-ladderize(pbtree(n=10, scale=10))
write.tree(tree, file="newick2.tre")