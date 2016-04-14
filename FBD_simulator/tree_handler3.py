import dendropy
import dendropy as dp
import numpy as np
from numpy import *
import scipy.stats
import time, csv
np.set_printoptions(suppress=True) # prints floats, no scientific notation
np.set_printoptions(precision=5)   # rounds all array elements to 3rd digit

## GLOBAL STUFF
rho = 1.0 # sampling
n_iterations = 10000
print_freq = 1000
sample_freq = 10
print_trees = False
runRJMCMC = False
tree_file = "/Users/daniele/Dropbox-personal/Dropbox/fossilizedBM/fossilBM/FBD_simulator/newick3.tre"

out_file_name="FBD.log" # % (dataset,args.j,model_name,clade_name,beta_value)
logfile = open(out_file_name , "wb") 
wlog=csv.writer(logfile, delimiter='\t')
head="it\tposterior\tlikelihood\tprior\tL\tM\tQ\tI"
wlog.writerow(head.split('\t'))
logfile.flush()

###### PRIORS
def prior_gamma(L,a=1,b=1): 
	return scipy.stats.gamma.logpdf(L, a, scale=1./b,loc=0)
def prior_normal(L,sd): 
	return scipy.stats.norm.logpdf(L,loc=0,scale=sd)
def prior_cauchy(x,s):
	return scipy.stats.cauchy.logpdf(x,scale=s,loc=0)


###### READ DATA 
tree = dp.Tree.get(path=tree_file, schema="newick")

# find fossils
tip_dist_from_root,taxa_labels=list(),list()
nd = tree.leaf_node_iter()
for i, tip in enumerate(nd):
	#print tip.distance_from_root(), tip.taxon._label
	tip_dist_from_root.append(tip.distance_from_root())
	taxa_labels.append(tip.taxon._label)


tip_dist_from_root=np.array(tip_dist_from_root)
taxa_labels = np.array(taxa_labels)
tree_root_age = max(tip_dist_from_root)
fossil_labels=taxa_labels[tip_dist_from_root<(tree_root_age-0.00001)]
m_fossils = len(fossil_labels)  
print m_fossils, "fossils:",fossil_labels


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

	#lik1 = 1./( lam*(1-p0hat(x1)) )**2
	#lik2 = 4*lam*rho / q(x1)
	#lik3 = np.prod( 4*lam*rho/q(x) )
	#lik4 = np.prod( psi*g * ( 2*lam*p0(y)*q(y)/q(z) )**I )

	log_lik1 = -(log(lam)+log(1-p0hat(x1)) )*2
	log_lik2 = log(4*lam*rho) - log(q(x1))
	log_lik3 = np.sum( log(4*lam*rho)-log(q(x)) )
	log_lik4 = np.sum( log(psi*g) + ( log(2*lam) + log(p0(y))+ log(q(y))-log(q(z)) )*I )
	#print sum(log(lik1)+log(lik2)+log(lik3)+log(lik4))- (log_lik1+log_lik2+log_lik3+log_lik4)

	return (log_lik1+log_lik2+log_lik3+log_lik4)

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
	if node.parent_node is None: return False
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


def move_fossil_tip_RJ(tip,set_mrca=False):
	# MRCA= tree.find_node_with_taxon_label("t2")
	# get distance from root of stem of fossil 
	br_length = tip.edge_length
	tip_age =tip.distance_from_root()
	root_dist_fossil_stem  = tip.distance_from_root()-br_length
	#print "NODE1:",tip._get_parent_node()
	# get root dist fossil tip
	fossil_tip_dist_root = tip.distance_from_root()
	remove_fossil(tip)
	# get mrca
	MRCA = set_mrca		
	root_dist_mrca= MRCA.distance_from_root()
	stem_mrca= MRCA._get_parent_node()
	root_dist_stem_mrca = stem_mrca.distance_from_root()
	len_stem_mrca= stem_mrca.edge_length
	#tip.edge_length=0.0000000001
	add_fossil(MRCA, tip)
	anc=tip._get_parent_node()
	root_dist_temp_node=anc.distance_from_root()
	#print dist_from_new_stem, root_dist_mrca-tip_age, root_dist_mrca,root_dist_stem_mrca
	alter= root_dist_temp_node-root_dist_fossil_stem  #+(10-root_dist_mrca)
	alter_node_age(tip, alter)




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
				if min(abs(rec_node_ages-age))<0.00001: 
					#node.label=""
					ind_x_z.append(0)
				else:
					#node.label="fossil"
					ind_x_z.append(1)
					fossil_node_list.append(node)
			else:
				if node.label=="fossil": 
					is_fossil=1
					fossil_node_list.append(node)
				ind_x_z.append(is_fossil)
				
			dists.append(age)
			
	ages=np.array(dists)
	#tree.print_plot(plot_metric='length',show_internal_node_labels=1) 
	#print "indexes", ind_x_z
	return ages,np.array(ind_x_z),fossil_node_list


def calc_gamma(tip, set_br_length=False):
	if set_br_length is False: br_length = tip.edge_length
	else: br_length=set_br_length
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
	
	if set_br_length is False:
		return len(candidate_nodes)
	else: # use in the delete move (RJMCMC)
		return candidate_nodes



####### START MCMC	
for it in range(n_iterations):
	hasting  = 0
	jacobian = 0
	
	if it==0:
		### MAKE a few TIPS BECOME FOSSIL
		#fossil_labels = ["t1","t2","t4","t9"] #["t2"] #
		fossil_tip_list = list()
		for f in fossil_labels: 
			fossil_tip=tree.find_node_with_taxon_label(f)
			fossil_edge = fossil_tip.edge
		#	fossil_edge.length = fossil_edge.length*0.15
			fossil_tip_list.append(fossil_tip)
			parent = fossil_tip._get_parent_node()
			parent.label="fossil"
			print f, fossil_tip
		
		if print_trees is True:
			print "Original tree:"
			tree.print_plot(plot_metric='length',show_internal_node_labels=1)
		
		# get node ages and identifiers for fossil nodes vs extant nodes
		n_ages, index_node,fossil_int_nodes = mod_calc_node_ages_(tree)
		reconstructed_node_ages = n_ages[index_node==0]# these node ages won't be touched
		# checks
		# nd= tree.postorder_node_iter(filter_fn=exclude_root) #  
		# original_root_dists= np.array([n.distance_from_root() for n in nd])
		# n1=  tree.mrca(taxon_labels=["t5","t3"])
		# orig_n1_age= n1.distance_from_root()
		# a,b,fossil_int_nodes =mod_calc_node_ages_(tree,rec_node_ages=reconstructed_node_ages)
		#
		# get ages of fossil tips
		tip_ages=np.zeros(len(fossil_labels))
		for i, f in enumerate(fossil_labels):
			tip=tree.find_node_with_taxon_label(f)
			tip_ages[i]=tree_root_age-tip.distance_from_root()
		
		par=[0.2,0.1,0.1]
		parA=par
		# indicators ancestral fossils = 0, fossil tips = 1
		Ind = np.ones(len(fossil_labels))
		
		
	rr=np.random.random(3)
	#if it<3: print it, rr, fossil_int_nodes
	# MOVE FOSSILS (ADD CONSTRAINTS!)
	if rr[0]<.4:
		# update attachment time
		r_fossil = np.random.randint(0,len(fossil_labels))
		#print "Before tree (%s):" % (fossil_labels[r_fossil])
		#tip_mod=fossil_tip_list[r_fossil]
		tip_mod=tree.find_node_with_taxon_label(fossil_labels[r_fossil])
		# update attachment placement
		alter_node_age(tip_mod)
		move_fossil_tip(tip_mod,3)
		#print gamma_f
	# UPDATE FOSSIL ATTACHMENT AGE (ADD CHECK FOR 0 BR LENGTH)
	elif rr[0]<0.8:
		# only update branches with non-zero br lengths
		# in theory internal nodes should all have >0 br length, because 
		# only fossil-tips can be assigned to a branch (hence having a 0 br length) 
		r_fossil_node = np.random.randint(0,len(fossil_int_nodes))
		alter_node_age(fossil_int_nodes[r_fossil_node],int_node=True)
	elif rr[0]<0.95 and runRJMCMC is True:
		# m_fossils: number of fossils
		# k_fossils: number of fossils that are ancestral
		index_k_fossils = (Ind==0).nonzero()[0]
		k_fossils = len(index_k_fossils)
		# ADD BRANCH with PROB g
		if k_fossils < m_fossils: g = 0.5 
		if m_fossils==k_fossils:  g = 1   # all are ancestral
		elif k_fossils==0:        g = 0   # none ancestral
		if g>rr[1]:
			r_index = np.random.choice(index_k_fossils)
			r_fossil_taxon = tree.find_node_with_taxon_label(fossil_labels[r_index])
			y_f = tip_ages[r_index] # age of the fossil
			anc_fossil_1 = r_fossil_taxon.parent_node # ancestor of the zero br length node
			# x_i = anc_fossil_1.parent_node # ancestor of the fossil
			u = rr[2]
			x_i_y_f = 0+anc_fossil_1.edge_length
			new_br_length = x_i_y_f*u # sample random br length
			alter_node_age(anc_fossil_1,alter=new_br_length,int_node=True)
			new_Ind = np.zeros(len(fossil_labels))+Ind
			new_Ind[r_index] = 1 # anc fossil is now a tip fossil
			new_index_k_fossils = (Ind==0).nonzero()[0]
			new_k_fossils = len(new_index_k_fossils)
			# hasting ratio
			if m_fossils==k_fossils and new_k_fossils>0: hasting = 0.5
			elif new_k_fossils==2:                       hasting = 2.0
			else:                                        hasting = 1.0
			# Jacobian
			jacobian = x_i_y_f
		# DELETE BRANCH (set br length to 0)
		else:
			tree.print_plot(plot_metric='length',show_internal_node_labels=1)
			index_tip_fossils = (Ind==1).nonzero()[0]
			r_index = 1 # np.random.choice(index_tip_fossils)
			r_fossil_taxon = tree.find_node_with_taxon_label(fossil_labels[r_index])
			# get possible attachment
			cand_lineages = calc_gamma(r_fossil_taxon, set_br_length=0)
			r_lineage = cand_lineages[1] # cand_lineages[np.random.randint(0,len(cand_lineages))]
			
			alter_node_age(r_fossil_taxon)
			move_fossil_tip_RJ(r_fossil_taxon,set_mrca=r_lineage)
			
	else:
		par, hasting = update_multiplier_proposal(parA,1.1)
	
		
	
	
	
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
		
	prior = sum(prior_gamma(par))
	lik= calc_FBD_lik(x=node_ages,z=fossil_sp_time,y=tip_ages,g=gamma_f,I=Ind,rates=par)
	
	if it==0: 
		likA=lik
		priorA=prior
		postA=likA+priorA
		IndA=Ind
	
	if (lik + prior) - postA + hasting >= log(np.random.random()):
		postA=lik+prior
		likA=lik
		priorA=prior
		parA=par
		IndA=Ind
	

	if it % print_freq==0:
		print it, likA, parA
		if print_trees is True: tree.print_plot(plot_metric='length',show_internal_node_labels=1) #
	
	if it % sample_freq==0:
		log_state=[it,postA,likA,priorA]+list(parA)+[np.mean(IndA)]
		wlog.writerow(log_state)
		logfile.flush()
	#time.sleep(0.2)






#nd= tree.postorder_node_iter(filter_fn=exclude_root) #  
#print np.sort(original_root_dists)
#print np.sort(np.array([n.distance_from_root() for n in nd]))
#nd= tree.postorder_node_iter(filter_fn=exclude_root) #  
#print np.sort(original_root_dists) - np.sort(np.array([n.distance_from_root() for n in nd]))





