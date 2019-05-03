import random
# given a network size, generates a uniform load
def generate_uniform_traffic(nnodes):
	load = 1. / float(nnodes * (nnodes - 1))
	tm = [0.] * nnodes
	for i in range(nnodes):
		tm[i] = [load] * nnodes
		tm[i][i] = 0.	
	return tm

# generates a traffic matrix such that a nnode only communicates with nodes 
# that are within k unit indices of one-another
def generate_closest_neighbor_traffic(nnodes, k):
	tm = [0. ] * nnodes 
	for i in range(nnodes):
		tm[i] = [0.] * nnodes
	return tm


def generate_bipartite_traffic(nnodes):
	tm = [0. ] * nnodes 
	for i in range(nnodes):
		tm[i] = [0.] * nnodes
	return tm

# returns a tm where tm is split into regions of clusters which are highly 
# communicative, and no communication elsewhere
def generate_clustered_traffic(nnodes, nclusters):
	tm = [0. ] * nnodes 
	num_nodes_per_cluster = nnodes / nclusters
	for i in range(nnodes):
		tm[i] = [0.] * nnodes
	return tm

# generates the traffic crossing from one dimension of the hyper x plan to the next in the final dimension
# group to group is generated in a bipartite like manner
def generate_hpx_intergroup_bipartite(tapered_hyperx):
	num_hyperx_groups = tapered_hyperx.num_groups()
	num_intragroup_switches = tapered_hyperx.num_intragroup_switches()
	switches_by_groups = [0] * num_hyperx_groups
	for i in range(num_hyperx_groups):
		switches_by_groups[i] = []
	# insert all the switches into the set
	for coord in tapered_hyperx.adjacency_list.keys():
		switches_by_groups[coord[-1]].append(tuple(coord))
	groups = range(num_hyperx_groups)
	intergroup_tm = [0] * num_hyperx_groups
	# if there are even number of groups
	tm = [0] * (num_hyperx_groups * num_intragroup_switches)
	for swid in range(num_hyperx_groups * num_intragroup_switches):
		tm[swid] = [0] * (num_hyperx_groups * num_intragroup_switches)
	# now finally distribute the traffic across its constituents
	for src_group in groups:
		src_set = (src_group <= num_hyperx_groups/2)
		for dst_group in groups:
			dst_set = (dst_group <= num_hyperx_groups/2)
			if src_set != dst_set:
				for i in range(num_intragroup_switches):
					coord = switches_by_groups[src_group][i]
					target_coord = switches_by_groups[dst_group][random.randint(0, len(switches_by_groups[dst_group]) - 1)]
				 	src_id = tapered_hyperx.coordinates_to_id[coord]
				 	dst_id = tapered_hyperx.coordinates_to_id[target_coord]
				 	tm[src_id][dst_id] += 1000. / num_intragroup_switches
	return tm


	
