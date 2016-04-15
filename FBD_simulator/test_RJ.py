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


fossil_tip_list = list()
for f in fossil_labels: 
	fossil_tip=tree.find_node_with_taxon_label(f)
	fossil_edge = fossil_tip.edge
#	fossil_edge.length = fossil_edge.length*0.15
	fossil_tip_list.append(fossil_tip)
	parent = fossil_tip._get_parent_node()
	parent.label="fossil"
	print f, fossil_tip

print "Original tree:"
tree.print_plot(plot_metric='length',show_internal_node_labels=1)

treeA = tree.clone(depth=1)

############################################################################
#			CODE TO REMOVE A FOSSIL BRANCH
############################################################################
# get node ages and identifiers for fossil nodes vs extant nodes
n_ages, index_node,fossil_int_nodes = mod_calc_node_ages_(tree)
reconstructed_node_ages = n_ages[index_node==0]# these node ages won't be touched
# get ages of fossil tips
tip_ages=np.zeros(len(fossil_labels))
for i, f in enumerate(fossil_labels):
	tip=tree.find_node_with_taxon_label(f)
	tip_ages[i]=tree_root_age-tip.distance_from_root()

par=[0.2,0.1,0.1]
parA=par



index_tip_fossils = (Ind==1).nonzero()[0]
r_index = np.random.choice(index_tip_fossils)
r_fossil_taxon = tree.find_node_with_taxon_label(fossil_labels[r_index])
print "\n\n\n FOSSIL CHANGED:", r_fossil_taxon
cand_lineages = calc_gamma(r_fossil_taxon, set_br_length=0)
r_lineage = cand_lineages[np.random.randint(0,len(cand_lineages))]
r_lineage
r_fossil_taxon
move_fossil_tip_RJ(r_fossil_taxon,set_mrca=r_lineage)
tree.print_plot(plot_metric='length',show_internal_node_labels=1)
Ind[r_index]=0


############################################################################
#			CODE TO ADD A FOSSIL BRANCH
############################################################################
index_k_fossils = (Ind==0).nonzero()[0]
k_fossils = len(index_k_fossils)
# ADD BRANCH with PROB g
if k_fossils < m_fossils: 
	g = 0.5 

if m_fossils==k_fossils:  
	g = 1   # all are ancestral
elif k_fossils==0:        
	g = 0   # none ancestral

r_index = np.random.choice(index_k_fossils)
r_fossil_taxon = tree.find_node_with_taxon_label(fossil_labels[r_index])
print "\n\n\n FOSSIL CHANGED:", r_fossil_taxon
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

jacobian = x_i_y_f
tree.print_plot(plot_metric='length',show_internal_node_labels=1)
Ind[r_index]=1
