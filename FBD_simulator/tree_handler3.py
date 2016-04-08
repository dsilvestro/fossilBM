import dendropy
import dendropy as dp
import numpy as np
from numpy import *
import time
np.set_printoptions(suppress=True) # prints floats, no scientific notation
np.set_printoptions(precision=5)   # rounds all array elements to 3rd digit

## GLOBAL STUFF
tree_root_age = 10.
rho = 1.0 # sampling

################### START LIKELIHOOD FUNCTIONS ###################
def calc_FBD_lik(x,z,y,g,I,rates=[0.5,0.2,0.1]):
	# x = recon_node_ages
	x1=max(x)
	# z = fossil_attach_age
	# y = fossil_age
	# g = possible_attachments_per_fossil
	# I = fossil_indicator # 0: along-branch, 1: fossil tip
	lam = rates[0] # spec rate
	mu  = rates[1] # ex rate
	psi = rates[2] # preservation


	c1= np.abs( np.sqrt( (lam-mu-psi)**2 + 4*lam*psi) )
	c2 = -(lam-mu-2*lam*rho-psi)/c1
	
	def q(t):
		return 2 * (1-c2**2) + exp(-c1*t) * (1-c2)**2 + exp(c1*t) * (1+c2)**2 
	
	def p0(t):
		return 1 + ( -(lam-mu-psi) + c1 * (exp(-c1*t)*(1-c2)-(1+c2))/(exp(-c1*t)*(1-c2)+(1+c2)) ) / 2*lam

	def p0hat(t):
		return 1 - (rho*(lam-mu)) / ( lam*rho + (lam*(1-rho)-mu) * exp(-(lam-mu)*t) ) 

	lik1 = 1./( lam*(1-p0hat(x1)) )**2
	lik2 = 4*lam*rho / q(x1)
	lik3 = np.prod( 4*lam*rho/q(x) )
	lik4 = np.prod( psi*g * ( 2*lam*p0(y)*q(y)/q(z) )**I )
	return sum(log(lik1)+log(lik2)+log(lik3)+log(lik4))
	

###################  END LIKELIHOOD FUNCTIONS  ###################
def update_multiplier_proposal(i,d):
	S=shape(i)
	u = np.random.uniform(0,1,S)
	l = 2*log(d)
	m = exp(l*(u-.5))
 	ii = i * m
	return ii, sum(log(m))


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

def alter_node_age(node,alter=0,int_node=False):
	if int_node is False: # node.is_leaf():
		anc_fossil = node._get_parent_node()
	else: anc_fossil=node
	try: 
		inc_edges = anc_fossil.get_incident_edges()
		#print anc_fossil
		#print node, "alter", alter
	except: 
		#print anc_fossil
		#print node, "alter", alter
		quit()
	# br lengths
	BRs= [i.length for i in inc_edges]
	#print "BR lengths:",BRs
	# min total br length (fossil tip and stem of it ancestor)
	minBR = min(BRs[0],BRs[1])
	if alter==0: 
		alter = round(np.random.uniform(-BRs[2],minBR),3)
		# alter 2 desc and stem
		inc_edges[0].length -= alter
		inc_edges[1].length -= alter
		inc_edges[2].length += alter
	else:
		inc_edges[0].length += alter
		inc_edges[1].length += alter
		inc_edges[2].length -= alter		
	#print alter, -BRs[2],minBR

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
	#return len(candidate_nodes) # number of possible attachmnents (gamma in FBD likelihood)



def mod_calc_node_ages(tree, ultrametricity_precision=False, internal_only=False):
	"""
	MODIFIED FROM DENDROPY: v
	Adds an attribute called "age" to  each node, with the value equal to
	the sum of edge lengths from the node to the tips. If the lengths of
	different paths to the node differ by more than ``ultrametricity_precision``, then a
	ValueError exception will be raised indicating deviation from
	ultrametricity. If ``ultrametricity_precision`` is negative or False, then this check
	will be skipped.
	"""
	ages = []
	for node in tree.postorder_node_iter():
		ch = node.child_nodes()
		if len(ch) == 0:
			node.age = 0.0
			if not internal_only:
				ages.append(node.age)
		else:
			first_child = ch[0]
			node.age = first_child.age + first_child.edge.length
			#if 2>1: #not (ultrametricity_precision is None or ultrametricity_precision is False or ultrametricity_precision < 0):
			#	for nnd in ch[1:]:
			#		try:
			#			ocnd = nnd.age + nnd.edge.length
			#		except TypeError:
			#			nnd.edge.length = 0.0
			#			ocnd = nnd.age
			#		d = abs(node.age - ocnd)
			#		if  d > ultrametricity_precision:
			#			# raise ValueError("Tree is not ultrametric. Node '{}': expecting {}, but found {}".format(node.label, node.age, ocnd))
			#			raise error.UltrametricityError("Tree is not ultrametric within threshold of {threshold}: {deviance}".format(
			#				threshold=ultrametricity_precision,
			#				deviance=d,
			#				))
			ages.append(node.age)
	return ages



def mod_calc_node_ages_(tree,root_age=tree_root_age,rec_node_ages=[]):
	"""
	Adds attribute "root_distance" to each node, with value set to the
	sum of edge lengths from the node to the root. Returns list of
	distances. If ``return_leaf_distances_only`` is True, then only
	leaf distances will be true.
	"""
	dists = []
	ind_x_z = []
	fossil_node_list=[]
	for node in tree.preorder_node_iter():
		is_fossil=0
		if node._parent_node is None:
			node.root_distance = 0.0
		else:
			node.root_distance = node.edge.length + node._parent_node.root_distance
		if node.is_leaf() is False:
			age = root_age-node.root_distance
			if len(rec_node_ages)>0:
				if min(abs(rec_node_ages-age))<0.00000001: 
					#node.label=""
					ind_x_z.append(0)
				else:
					#node.label="fossil"
					ind_x_z.append(1)
					fossil_node_list.append(node)
			else:
				if node.label=="fossil": is_fossil=1
				ind_x_z.append(is_fossil)
				
			dists.append(age)
			
	ages=np.array(dists)
	#tree.print_plot(plot_metric='length',show_internal_node_labels=1) 
	#print "indexes", ind_x_z
	return ages,np.array(ind_x_z),fossil_node_list


def calc_gamma(tip):
	br_length = tip.edge_length
	tip_age =tip.distance_from_root()
	root_dist_fossil_stem  = tip.distance_from_root()-br_length
	# get root dist fossil tip
	fossil_tip_dist_root = tip.distance_from_root()
	#print root_dist_fossil_stem, fossil_tip_dist_root
	nd= tree.postorder_node_iter(filter_fn=exclude_root) #  
	candidate_nodes=[]
	for n, node in enumerate(nd): 
		# get parent
		anc_fossil = node.parent_node
		t1 = anc_fossil.distance_from_root()
		t2 = node.distance_from_root()
		#print n, t1,t2,
		diff = np.array([(t1-root_dist_fossil_stem), t2-fossil_tip_dist_root]) # a fossil cannot be assigned to itself
		if t1<root_dist_fossil_stem and t2>root_dist_fossil_stem and max(abs(diff))>0.00001:
			candidate_nodes.append(node)
			#print "*",diff 
		else: pass #print "",diff
	return len(candidate_nodes)


###### READ DATA 
tree = dp.Tree.get(path="/Users/daniele/Dropbox-personal/Dropbox/fossilizedBM/FBD_simulator/newick2.tre", schema="newick")
	
for it in range(10000):
	#tree =  dp.Tree.get(data="[&R] ((((((t7:0.3180281068,t8:0.3180281068):3.627021941,t4:3.945050048):0.08067589546,t3:4.025725943):2.958874374,t2:6.984600317):2.129306342,((t5:0.6972928451,t6:0.6972928451):3.308282054,(t9:0.06071668103,t10:0.06071668103):3.944858218):5.10833176):0.8860933409,t1:10);", schema="newick")
	hasting=0
	if it==0:
		### MAKE a few TIPS BECOME FOSSIL
		fossil_labels = ["t1","t2","t4","t9"] #["t2"] #
		fossil_tip_list = list()
		for f in fossil_labels: 
			fossil_tip=tree.find_node_with_taxon_label(f)
			fossil_edge = fossil_tip.edge
			fossil_edge.length = fossil_edge.length*0.15
			fossil_tip_list.append(fossil_tip)
			parent = fossil_tip._get_parent_node()
			parent.label="fossil"
		print "Original tree:"
		tree.print_plot(plot_metric='length',show_internal_node_labels=1)
		# get node ages and identifiers for fossil nodes vs extant nodes
		n_ages, index_node,fossil_int_nodes = mod_calc_node_ages_(tree)
		reconstructed_node_ages = n_ages[index_node==0]# these node ages won't be touched
		# checks
		nd= tree.postorder_node_iter(filter_fn=exclude_root) #  
		original_root_dists= np.array([n.distance_from_root() for n in nd])
		n1=  tree.mrca(taxon_labels=["t5","t3"])
		orig_n1_age= n1.distance_from_root()
		a,b,fossil_int_nodes =mod_calc_node_ages_(tree,rec_node_ages=reconstructed_node_ages)
		# get ages of fossil tips
		tip_ages=np.zeros(len(fossil_labels))
		for i, f in enumerate(fossil_labels):
			tip=tree.find_node_with_taxon_label(f)
			tip_ages[i]=tree_root_age-tip.distance_from_root()
		par=[0.5,0.2,0.1]
		parA=par
		
	rr=np.random.random()
	if rr<.5:
		# update attachment time
		r_fossil = np.random.randint(0,len(fossil_labels))
		#print "Before tree (%s):" % (fossil_labels[r_fossil])
		#tip_mod=fossil_tip_list[r_fossil]
		tip_mod=tree.find_node_with_taxon_label(fossil_labels[r_fossil])
		# update attachment placement
		alter_node_age(tip_mod)
		move_fossil_tip(tip_mod,3)
		#print gamma_f
	elif rr<0.95:
		r_fossil_node = np.random.randint(0,len(fossil_int_nodes))
		alter_node_age(fossil_int_nodes[r_fossil_node],int_node=True)
	else:
		par, hasting = update_multiplier_proposal(parA,1.2)
		
	# get vector of speciation times
	# rename node
	#tree.print_plot(plot_metric='length',show_internal_node_labels=1) #
	#tip_mod=tree.find_node_with_taxon_label(fossil_labels[r_fossil])
	#parent = tip_mod._get_parent_node()
	#parent.label="fossil"
	#print "After tree"
	
	#print reconstructed_node_ages
	a,b,fossil_int_nodes =mod_calc_node_ages_(tree,rec_node_ages=reconstructed_node_ages)
	node_ages = a[b==0]
	fossil_sp_time = a[b==1]
	
	
	
	gamma_f = np.zeros(len(fossil_labels))
	for i, f in enumerate(fossil_labels):
		tip=tree.find_node_with_taxon_label(f)
		gamma_f[i] = calc_gamma(tip)
	
	lik= calc_FBD_lik(x=node_ages,z=fossil_sp_time,y=tip_ages,g=gamma_f,I=np.ones(len(fossil_labels)),rates=par)
	if it==0: 
		likA=lik
		prior=0
		postA=likA+prior
	
	if (lik + prior) - postA + hasting >= log(np.random.random()):
		postA=lik+prior
		likA=lik
		priorA=prior
		parA=par
	

	if it % 1000==0:
		print it, likA, parA
		#tree.print_plot(plot_metric='length',show_internal_node_labels=1) #
	
	#time.sleep(0.2)






#nd= tree.postorder_node_iter(filter_fn=exclude_root) #  
#print np.sort(original_root_dists)
#print np.sort(np.array([n.distance_from_root() for n in nd]))
#nd= tree.postorder_node_iter(filter_fn=exclude_root) #  
#print np.sort(original_root_dists) - np.sort(np.array([n.distance_from_root() for n in nd]))


n1=  tree.mrca(taxon_labels=["t5","t3"])# tree.find_node_with_taxon_label("t5"), tree.find_node_with_taxon_label("t3"))
print n1.distance_from_root()




