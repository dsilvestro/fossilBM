import sys
from time import time
import dendropy as dp
from numpy.random import poisson, uniform

BD_treefile=sys.argv[1]
FBD_treefile=sys.argv[2]
q=float(sys.argv[3])

# commandline: python simulate_fossils_on_tree.py BD_treefile FBD_treefile q

#####################################################################################################
# This script takes a tree with extant and extinct lineages. 
# It adds poisson numbers of fossil on each branch, which depend on the branch length and the rate 
# of fossil discovery, q. Each fossil takes the form of a new leaf having a 0.00001 branch length
# In a second step, extinct lineages (except fossils) are pruned from the tree. The resulting tree 
# can therefore have a smaller root age than the original tree. In the case where the root node of 
# the tree has a branch length, fossils older than root can be added. 
#####################################################################################################

def add_fossil_tip(root_age, split_node, fossil_tip, div_age):
	""" Function which add a tip in the branch upstream of a given node (split node).
	1. The split node is cloned (clone_node). The cloned node encompasses therefore the same children than the split node.
	2. The idea is that the combined edge lengths of the clone node and split node should be equal 
		to the original split node edge length. Therefore, the clone node edge length corresponds to original edge length minus
		the divergence age (div_age) of the fossil. And we substract the clone node edge length to the original split node.
	3. Then we kill the children of the split node (all the substree linked to it) and set its taxon object to None (important if it was a leaf).
	4. Finally, we add the clone node (with the subtree linked to it) AND a new tip (fossil tip) to the the childless split node."""
	# get the clone of the split node
	clone_node=split_node.clone(depth=2)
	# edge length of the clone node
	clone_node.edge_length=split_node.distance_from_root()-(root_age-div_age)
	# substract edge length of the clone node to the original split node
	split_node.edge_length-=clone_node.edge_length
	# kill its childs
	split_node.clear_child_nodes()
	# remove taxon
	split_node.taxon=None
	# add the fossil_tip and the clone node to the childless ( :-( ) split node
	split_node.add_child(clone_node)
	split_node.add_child(fossil_tip)

#####################################################################################################

def is_extant_leaf(tree, node):
	""" Check if a leaf is extant and not extinct """
	if node.is_leaf() and node.distance_from_root()>(tree.max_distance_from_root()+tree.seed_node.edge_length-0.00001):
		return True
	else: return False

#####################################################################################################

def simulate_fossils_on_tree(tree,q):
	""" the main function to simulate fossils on a tree """
	# update its bipartitions (way of indexing the edges)
	tree.update_bipartitions()
	# get age of the root
	root_age=tree.max_distance_from_root()
	# store taxon label needed to be kept (the extant and fossils leaves)
	taxon_label_to_keep=[]
	# fossil id
	f=1
	# loop through bipartitions
	for i in range(len(tree.bipartition_encoding)):
		# access the edge through the bipartitions-edges hash
		edge=tree.bipartition_edge_map[tree.bipartition_encoding[i]]
		# we add the tip node, if it is an extant leaf, to the taxon_label_to_keep
		if is_extant_leaf(tree,edge.head_node)==True:
			taxon_label_to_keep.append(edge.head_node.taxon.label)
		# generate a poisson number of fossil 
		F=poisson(q*edge.length)
		# get F divergence age values. 
		# They need to be sorted in order to use always the same split node in the function 'add_fossil_tip'.
		F_div_ages=sorted([uniform(0,edge.length)+(root_age-edge.head_node.distance_from_root()) for x in range(F)])
		# loop through fossil divergence ages
		for i in range(F):
			# create a new fossil tip object, with a very small branch length
			fossil_tip=dp.Node(edge_length=0.00001, taxon=dp.Taxon(label='f%s'%(f)))
			# add the fossil
			add_fossil_tip(root_age,edge.head_node,fossil_tip,F_div_ages[i])
			# add the fossil tip taxon label to the taxon_label_to_keep
			taxon_label_to_keep.append('f%s'%(f))
			f+=1	
	# write temp_tree with fossil and extinct lineages
	# tree.write(path="10/%s_temp.tree"%(str(iterator)), schema="newick")
	# update taxonnamespace
	tree.update_taxon_namespace()
	# Now, we prune extinct tips
	tree.retain_taxa_with_labels(taxon_label_to_keep)
	return tree

### main function ##################################################################################

def main_IO_BDtree_2_FBDtree(BD_treefile,FBD_treefile,q):
	""" main function that takes a tree file and write the FBD tree """
	# read the tree (comming from treeSim)
	tree = dp.Tree.get(path=BD_treefile, schema="newick",rooting='force-rooted')
	# add fossils
	tree = simulate_fossils_on_tree(tree,q)
	# write fbd tree
	tree.write(path=FBD_treefile, schema='newick')

#####################################################################################################

main_IO_BDtree_2_FBDtree(BD_treefile,FBD_treefile,q)

# Benchmarking :
# import glob
#paths=glob.glob('./1000/*.nwk')
##paths=['tree.nwk']
#print('start!')
#for i in range(len(paths)):
#	print str(i)
#	start=time()
#	# read tree
#	tree = dp.Tree.get(path=paths[i], schema="newick",rooting='force-rooted')
#	# number edges
#	n_edges=len(tree.edges())
#	tree = simulate_fossils_on_tree(tree,q)
#	end=time()-start
#	print("%s\n%s\tedges\n%s"%(paths[i],str(n_edges),str(end)))
#	# write final tree
#	tree.write(path="1000/%s_out.tree"%(str(i)), schema="newick")
#print('end!')


