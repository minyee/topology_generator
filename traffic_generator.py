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
def generate_hpx_intergroup_clustered(tapered_hyperx, nclusters):
	num_groups = tapered_hyperx.S[-1]
	num_groups_per_cluster = num_groups / nclusters
	leftover_groups = num_groups % nclusters
	groups_in_cluster = [[]] * nclusters
	group_id = 0 
	intergroup_tm = [0] * num_groups
	for group in range(num_groups):
		intergroup_tm[group] = [0] * num_groups
	for cluster in range(nclusters):
		for i in range(num_groups_per_cluster):
			groups_in_cluster[cluster].append(group_id)
			group_id += 1
		if leftover_groups > 0:
			leftover_groups -= 1
			groups_in_cluster[cluster].append(group_id)
			group_id += 1
	for cluster in range(nclusters):
		for src_group in groups_in_cluster[cluster]:
			for dst_group in groups_in_cluster[cluster]:
				if src_group != dst_group:
					intergroup_tm[src_group][dst_group] = 1.
	switches_in_group = {}
	tm = [0] * len(tapered_hyperx.adjacency_list.keys())
	for i in range(len(tapered_hyperx.adjacency_list.keys())):
		tm[i] = [0] * len(tapered_hyperx.adjacency_list.keys())
	for coord in tapered_hyperx.adjacency_list.keys():
		group = coord[-1]
		if group not in switches_in_group.keys():
			switches_in_group[group] = []
		switches_in_group[group].append(tapered_hyperx.coordinates_to_id[coord])
	#for cluster in nclusters:
	#	for src_group in groups_in_cluster[cluster]:
	#		for dst_group in groups_in_cluster[cluster]:
	#			if src_group != dst_group:
	#				intergroup_tm[src_group][dst_group] = tapered_hyperx.num_intragroup_switches() * 1000./len(groups_in_cluster[cluster] - 1)
	# finally, fill in the traffic matrix for switch to switch
	for cluster in range(nclusters):
		for src_group in groups_in_cluster[cluster]:
			for dst_group in groups_in_cluster[cluster]:
				if src_group != dst_group:
					for index in range(tapered_hyperx.num_intragroup_switches()):
						tm[switches_in_group[src_group][index]][switches_in_group[dst_group][index]] = 1000./ tapered_hyperx.num_intragroup_switches()
	return tm, intergroup_tm

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
	for i in range(num_hyperx_groups):
		intergroup_tm[i] = [0] * num_hyperx_groups
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
				intergroup_tm[src_group][dst_group] = 1.
				for i in range(num_intragroup_switches):
					coord = switches_by_groups[src_group][i]
					target_coord = switches_by_groups[dst_group][random.randint(0, len(switches_by_groups[dst_group]) - 1)]
				 	src_id = tapered_hyperx.coordinates_to_id[coord]
				 	dst_id = tapered_hyperx.coordinates_to_id[target_coord]
				 	tm[src_id][dst_id] += 1000. / num_intragroup_switches
	return tm, intergroup_tm

def read_traffic_trace_files(filename):
	tm_dict = {}
	tm = []
	with open(filename) as f:
		num_blocks = 0
		for line in f:
			values = line.split(' ')
			src = int(values[2])
			dst = int(values[3])
			size = int(values[4])
			num_blocks = max(num_blocks, src, dst)
			if src not in tm_dict:
				tm_dict[src] = {}
			if dst not in tm_dict[src]:
				tm_dict[src][dst] = size
			else:
				tm_dict[src][dst] += size
		num_blocks += 1
		num_blocks = max(len(tm_dict.keys()), num_blocks)
		#num_blocks = max(tm_dict.keys())
		tm = [0] * num_blocks
		for src in range(num_blocks):
			tm[src] = [0] * num_blocks
			if src not in tm_dict.keys():
				continue
			for dst in range(num_blocks):
				if dst not in tm_dict[src].keys():
					continue
				else:
					tm[src][dst] += tm_dict[src][dst]
	return tm

## reshapes the TM into a new_size tm
def reshape_tm(tm, new_size):
	old_size = len(tm)
	if old_size == new_size:
		return tm
	elif old_size > new_size:
		size_per_new_blocksize = old_size / new_size
		leftovers = old_size % new_size
		offset = 0
		new_blocks_includes = [0] * new_size # contains what the old block ids belonging to each new block
		old_blocks_belong = {}
		# figure out which old blocks belong to which new blocks
		for i in range(new_size):
			new_blocks_includes[i] = []
			for j in range(size_per_new_blocksize):
				new_blocks_includes[i].append(offset)
				old_blocks_belong[offset] = i
				offset += 1
			if leftovers > 0:
				leftovers -= 1
				new_blocks_includes[i].append(offset)
				old_blocks_belong[offset] = i
				offset += 1
		new_tm = [0] * new_size
		for src_new_block in range(new_size):
			new_tm[src_new_block] = [0] * new_size
		for i in range(old_size):
			new_block_src = old_blocks_belong[i]
			for j in range(old_size):
				if i != j:
					new_block_dst = old_blocks_belong[j]
					new_tm[new_block_src][new_block_dst] += tm[i][j]
		return new_tm
	else: 
		print "Unimplemented"
		return None
