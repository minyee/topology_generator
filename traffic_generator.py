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


def generate_hpx_