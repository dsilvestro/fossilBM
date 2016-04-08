import dendropy
import dendropy as dp
import numpy as np

	
###### TRANSFORM TREE NODES
# parse newick
tree =  dp.Tree.get(data="[&R] ((((((t7:0.3180281068,t8:0.3180281068):3.627021941,t4:3.945050048):0.08067589546,t3:4.025725943):2.958874374,t2:6.984600317):2.129306342,((t5:0.6972928451,t6:0.6972928451):3.308282054,(t9:0.06071668103,t10:0.06071668103):3.944858218):5.10833176):0.8860933409,t1:10);", schema="newick")
#tree = dp.Tree.get(path="/Users/daniele/Dropbox-personal/Dropbox/fossilizedBM/FBD_simulator/newick.tre", schema="newick")
tree.print_plot(plot_metric='length')
#tree.is_rooted

### MAKE a TIP BECOME A FOSSIL
fossil_tip=tree.find_node_with_taxon_label("t2")
# distance from root
fossil_tip_dist_root = fossil_tip.distance_from_root()
# edge
fossil_edge = fossil_tip.edge
fossil_edge.length = fossil_edge.length*0.15
tree.print_plot(plot_metric='length')
# get anc node of fossil
anc_fossil = fossil_tip._get_parent_node()

# UPDATE age fossil anc
for i in range(0):
	# get incident edges
	inc_edges = anc_fossil.get_incident_edges()
	# br lengths
	BRs= [i.length for i in inc_edges]
	# min total br length (fossil tip and stem of it ancestor)
	minBR = min(BRs[0],BRs[1])
	alter = np.random.uniform(-BRs[2],minBR)
	# alter 2 desc and stem
	inc_edges[0].length -= alter
	inc_edges[1].length -= alter
	inc_edges[2].length += alter
	tree.print_plot(plot_metric='length')





def alter_node_age(node):
	anc_fossil = node._get_parent_node()
	inc_edges = anc_fossil.get_incident_edges()
	# br lengths
	BRs= [i.length for i in inc_edges]
	# min total br length (fossil tip and stem of it ancestor)
	minBR = min(BRs[0],BRs[1])
	alter = np.random.uniform(-BRs[2],minBR)
	# alter 2 desc and stem
	inc_edges[0].length -= alter
	inc_edges[1].length -= alter
	inc_edges[2].length += alter
	

tree =  dp.Tree.get(data="[&R] ((((((t7:0.3180281068,t8:0.3180281068):3.627021941,t4:3.945050048):0.08067589546,t3:4.025725943):2.958874374,t2:6.984600317):2.129306342,((t5:0.6972928451,t6:0.6972928451):3.308282054,(t9:0.06071668103,t10:0.06071668103):3.944858218):5.10833176):0.8860933409,t1:10);", schema="newick")
#tree = dp.Tree.get(path="/Users/daniele/Dropbox-personal/Dropbox/fossilizedBM/FBD_simulator/newick.tre", schema="newick")
tree.print_plot(plot_metric='length')


MRCA=tree.mrca(taxon_labels=["t2","t7"])

ch3 = dendropy.Node(edge_length=1)

# MRCA is a leaf
MRCA= tree.find_node_with_taxon_label("t2")
MRCA1=MRCA.clone(depth=2)
MRCA1.edge_length=ch3.edge_length
MRCA.add_child(ch3)
MRCA.add_child(MRCA1)
MRCA.edge_length-=ch3.edge_length





MRCA.insert_child(index = 0,node =ch3) # using index=2 new node goes along the stem of MRCA
tree.resolve_polytomies(limit=2, update_bipartitions=False, rng=None)
tree.print_plot(plot_metric='length')

alter_node_age(ch3)
tree.print_plot(plot_metric='length')







### OTHER STUFF
MRCA=tree.mrca(taxon_labels=["t5","t9"])
tip=tree.find_node_with_taxon_label("t8")
MRCA.insert_new_child(tip)
#MRCA.insert_child(0,fossil_tip )
tree.print_plot(plot_metric='length')


# get node ages
tree.calc_node_ages(ultrametricity_precision=False) # allow for non-ultram trees

tree.find_node_with_label("t2")


#tree.resolve_polytomies(update_bipartitions=True)
tree.calc_node_ages(ultrametricity_precision=0.001) # 
nd= tree.ageorder_node_iter(include_leaves=False, filter_fn=None, descending=True)
ages=list()
edges=list()
for n, node in enumerate(nd): 
	"loop nodes"



#tree.resolve_polytomies(update_bipartitions=True)
tree.calc_node_ages(ultrametricity_precision=0.001) # 
nd= tree.ageorder_node_iter(include_leaves=False, filter_fn=None, descending=True)
ages=list()
edges=list()
for n, node in enumerate(nd): 
	ages.append(node.age)
	edges.append(node.edge_length)
	print n, node



tree = dendropy.Tree.get(path="/Users/daniele/Dropbox-personal/Dropbox/fossilizedBM/FBD_simulator/newick.tre", schema="newick")
tree.print_plot(plot_metric='age')
N = tree.nodes()
N[1].edge_length = 8
tree.print_plot(plot_metric='age')
N[1].num_child_nodes()es.append(node.age)
	edges.append(node.edge_length)
	print n, node



tree = dendropy.Tree.get(path="/Users/daniele/Dropbox-personal/Dropbox/fossilizedBM/FBD_simulator/newick.tre", schema="newick")
tree.print_plot(plot_metric='length')
N = tree.nodes()
N[1].edge_length = 8
tree.print_plot(plot_metric='length')
N[1].num_child_nodes()


# change edge lengths
def transform_func(x): return x*np.random.uniform(0.1,2)

for n, edge in enumerate(tree.preorder_edge_iter()): 
	print n, edge.length
	if n>0:
		edge.length = transform_func(edge.length)

tree.print_plot(plot_metric='length')


# loop over internal nodes
for n, internal_node in enumerate(tree.postorder_internal_node_iter()): 
	print n, internal_node.distance_from_root(),internal_node.taxon #internal_node._child_nodes
	if n>0: internal_node.distance_from_root= transform_func(internal_node.distance_from_root())

tree.print_plot(plot_metric='length')


# get list of all attributes
dir(internal_node)


# loop over all nodes
for n, node in enumerate(tree.postorder_node_iter()): 
	str=node.taxon
	print n, node.distance_from_root(),str #internal_node._child_nodes
	# alter one tip
	print str, type(str)
	


