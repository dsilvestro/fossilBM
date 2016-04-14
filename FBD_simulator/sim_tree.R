setwd("~/Dropbox-personal/Dropbox/fossilizedBM/fossilBM/FBD_simulator")
library(phytools)

tree<-ladderize(pbtree(n=10, scale=10))
write.tree(tree, file="newick2.tre")



library(TreeSim)
t =sim.bd.taxa(n=20, 1, lambda=0.1, mu=0.025, frac = 1, complete = TRUE, stochsampling = FALSE)
plot(t[[1]])
write.tree(t[[1]], file="newick3.tre")
