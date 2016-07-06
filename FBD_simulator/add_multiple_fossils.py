import dendropy as dp
import random
#-------------------------------------------------------------------------------------------
# FUNCTION READ FOSSIL TABLE
#-------------------------------------------------------------------------------------------

def read_fossil_table(filename):
	""" input: table with column 1 as fossil name and column 2 as fossil age
	output: A list of fossil_names and a list of fossil_ages"""
	fossil_names=list()
	fossil_ages=list()
	with open(filename, 'r') as file:
		for line in file:
			fossil_names.append(line.split()[0].strip('\"'))
			fossil_ages.append(float(line.split()[1]))
	return fossil_names,fossil_ages

#-------------------------------------------------------------------------------------------
# FUNCTION GENERATE DIVERGENCE TIMES
#-------------------------------------------------------------------------------------------

def get_divergence_age(fossil_ages, tree):
	""" input1: list of fossil_ages
	input2: tree
	output: list of divergence age, which is a random age between the fossil age
	and root age """
	root_age=tree.max_distance_from_root()
	div_ages=list()
	for fossil_age in fossil_ages:
		div_ages.append(random.uniform(root_age, fossil_age))
	return div_ages

#-------------------------------------------------------------------------------------------
# FUNCTION GET CANDIDATE NODES
#-------------------------------------------------------------------------------------------

def get_candidates(tree, div_age, fossil_names):
	"""input1: tree
	input2: divergence age (of fossil)
	input3: list of fossil names
	The function returns a list of candidate nodes (the nodes downstream the div. age) to attach the fossil. 
	However, it doesn't return candidates, which are already added fossil tips.
	"""
	candidates=list()
	for n in tree:
		if n.is_leaf():
			if n.taxon.label not in fossil_names:
				if div_age>n.distance_from_tip() and (div_age-n.distance_from_tip())<n.edge_length:
					candidates.append(n)
		elif div_age>n.distance_from_tip() and (div_age-n.distance_from_tip())<n.edge_length:
			candidates.append(n)			
	return candidates

#-------------------------------------------------------------------------------------------
# FUNCTION ADD_FOSSIL_TIP
#-------------------------------------------------------------------------------------------

def add_fossil_tip(split_node, fossil_tip, div_age):
	""" input1: one node of the tree
	input2: """
	# get the clone of the split node
	clone_node=split_node.clone(depth=2)
	# edge length of the clone node
	clone_node.edge_length=div_age-split_node.distance_from_tip()
	# substract edge length of the clone node to the original split node
	split_node.edge_length-=clone_node.edge_length
	# kill its childs
	split_node.clear_child_nodes()
	# remove taxon
	split_node.taxon=None
	# add the fossil_tip and the clone node to the childless ( :-( ) split node
	split_node.add_child(clone_node)
	split_node.add_child(fossil_tip)


#------------------------------------------------------------------------------------------
# Read tree
tree = dp.Tree.get(path='/Users/daniele/Dropbox-personal/Dropbox/fossilizedBM/code_FBD/fossilBM/FBD_simulator/tree.tre', schema="newick")
# Read fossil table
fossil_names, fossil_ages=read_fossil_table('fossils.txt')
# get divergence ages of all fossils
div_ages=get_divergence_age(fossil_ages, tree)

# add multiple fossils (test for first 5)
for i, fossil_name in enumerate(fossil_names):
	# get candidate nodes
	c=get_candidates(tree, div_ages[i], fossil_names)
	# choose a random candidate
	split_node=random.choice(c)
	# fossil's tip length
	tip_length=div_ages[i]-fossil_ages[i]
	# create fossil tip
	fossil_tip=dp.Node(edge_length=tip_length, taxon=dp.Taxon(label=fossil_name))
	# add the fossil
	add_fossil_tip(split_node, fossil_tip, div_ages[i])
	# see
	print(tree.as_string('newick'))
	tree.write(path=str(i)+"_output.tree", schema="newick")
#-------------------------------------------------------------------------------------------






