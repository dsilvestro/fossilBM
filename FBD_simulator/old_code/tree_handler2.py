import dendropy
import dendropy as dp
import numpy as np
import time

def add_fossil(node,tip):
	if node.is_leaf():
		#node= tree.find_node_with_taxon_label("t2")
		node1=node.clone(depth=2)
		node1.edge_length=tip.edge_length
		node.add_child(tip)
		node.add_child(node1)
		node.edge_length-=tip.edge_length
	else:
		##tip = tree.find_node_with_taxon_label(tip_label)
		node.insert_child(index = 2,node =tip) # using index=2 new node goes along the stem of MRCA
		try: tree.resolve_polytomies(limit=2, update_bipartitions=True, rng=None)
		except: pass


def exclude_root(node):
	if node.parent_node is None:
		return False
	else: return True

def remove_fossil(node):
	parent = node._get_parent_node()
	parent.remove_child(node, suppress_unifurcations=True)

def alter_node_age(node,alter=0):
	anc_fossil = node._get_parent_node()
	inc_edges = anc_fossil.get_incident_edges()
	# br lengths
	BRs= [i.length for i in inc_edges]
	# min total br length (fossil tip and stem of it ancestor)
	minBR = min(BRs[0],BRs[1])
	if alter==0: 
		alter = np.random.uniform(-BRs[2],minBR)
		# alter 2 desc and stem
		inc_edges[0].length -= alter
		inc_edges[1].length -= alter
		inc_edges[2].length += alter
	else:
		inc_edges[0].length += alter
		inc_edges[1].length += alter
		inc_edges[2].length -= alter	
	print alter, -BRs[2],minBR

def move_fossil_tip(tip,ind=0):
	# MRCA= tree.find_node_with_taxon_label("t2")
	# get distance from root of stem of fossil 
	br_length = tip.edge_length
	tip_age =tip.distance_from_root()
	root_dist_fossil_stem  = tip.distance_from_root()-br_length
	#print "NODE1:",tip._get_parent_node()
	# get root dist fossil tip
	fossil_tip_dist_root = tip.distance_from_root()
	remove_fossil(tip)
	#tree.print_plot(plot_metric='length')
	# select candidate node
	nd= tree.postorder_node_iter(filter_fn=exclude_root) #  
	candidate_nodes=[]
	for n, node in enumerate(nd): 
		# get parent
		anc_fossil = node.parent_node
		t1 = anc_fossil.distance_from_root()
		t2 = node.distance_from_root()
		#print n, t1,t2
		if t1<root_dist_fossil_stem and t2>root_dist_fossil_stem:
			candidate_nodes.append(node)
	
	# random candidate node
	ind= np.random.randint(0,len(candidate_nodes))
	#print "NODE2:",candidate_nodes[ind], ind
	MRCA = candidate_nodes[ind]
	root_dist_mrca= MRCA.distance_from_root()
	stem_mrca= MRCA._get_parent_node()
	root_dist_stem_mrca = stem_mrca.distance_from_root()
	len_stem_mrca= stem_mrca.edge_length
	add_fossil(MRCA, tip)
	
	anc=tip._get_parent_node()
	root_dist_temp_node=anc.distance_from_root()
	
	
	#print dist_from_new_stem, root_dist_mrca-tip_age, root_dist_mrca,root_dist_stem_mrca
	alter= root_dist_temp_node-root_dist_fossil_stem  #+(10-root_dist_mrca)
	alter_node_age(tip, alter)
	tip.edge_length = br_length
	#print "ages", tip_age, tip.distance_from_root(), "alter:",alter
	return len(candidate_nodes) # number of possible attachmnents (gamma in FBD likelihood)



tree = dp.Tree.get(path="/Users/daniele/Dropbox-personal/Dropbox/fossilizedBM/FBD_simulator/newick2.tre", schema="newick")
	
for i in range(100):
	# distance from root
	#fossil_tip_dist_root = fossil_tip.distance_from_root()
	#tree =  dp.Tree.get(data="[&R] ((((((t7:0.3180281068,t8:0.3180281068):3.627021941,t4:3.945050048):0.08067589546,t3:4.025725943):2.958874374,t2:6.984600317):2.129306342,((t5:0.6972928451,t6:0.6972928451):3.308282054,(t9:0.06071668103,t10:0.06071668103):3.944858218):5.10833176):0.8860933409,t1:10);", schema="newick")
	if i==0:
		### MAKE a TIP BECOME A FOSSIL
		#print "fossil age 1",fossil_tip_dist_root
		# edge
		fossil_tip=tree.find_node_with_taxon_label("t2")
		fossil_edge = fossil_tip.edge
		fossil_edge.length = fossil_edge.length*0.15
		print "Original tree:"
		tree.print_plot(plot_metric='length')
	print "New tree:"
	# update attachment time
	alter_node_age(fossil_tip)
	# update attachment placement
	gamma_f = move_fossil_tip(fossil_tip,3)
	print gamma_f
	tree.print_plot(plot_metric='length') #
	time.sleep(1)









