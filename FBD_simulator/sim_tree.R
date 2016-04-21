setwd("~/Dropbox-personal/Dropbox/fossilizedBM/fossilBM/FBD_simulator")
library(phytools)

tree<-ladderize(pbtree(n=10, scale=10))
write.tree(tree, file="newick2.tre")



library(TreeSim)
t =sim.bd.taxa(n=40, 1, lambda=0.2, mu=0.1, frac = 1, complete = TRUE, stochsampling = FALSE)
plot(ladderize(t[[1]]))
write.tree(t[[1]], file="newick40_0201.tre", tree.names="[&R] ")
