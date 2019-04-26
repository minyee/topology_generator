import sys
from functools import reduce
import random
import copy

## normally describes just a static network without optical circuit switches, direct network because every switch is supposed to be attached to a terminal

def cart_product(set1, set2):
	cart_prod = []
	for element1 in set1:
		for element2 in set2:
			cart_prod.append(element1 + element2)

	return cart_prod

# splits the universe set into 2 subsets, A and B, such that A union B equals universe.
# need to receive copy 
def splitset(universe1):
	universe = list(universe1)
	size_A = random.randint(1, len(universe)/2 + 1)
	A = []
	B = []
	for i in range(size_A):
		index = random.randint(0, len(universe) - 1)
		if index >= len(universe):
			print "universe : {}   ,   index : {}".format(len(universe), index)
		assert(index < len(universe))
		A.append(universe[index])
		del universe[index]
	B = universe
	return A, B

class HyperX:
	def __init__(self, L_arg, S_arg, K_arg, T_arg):
		self.L = L_arg
		self.S = S_arg ## Note that 
		assert(len(self.S) == 1 or len(self.S) == L_arg)
		if len(self.S) == 1:
			self.S = [self.S[0]] * self.L
		self.K = K_arg
		self.T = T_arg
		self.adjacency_list = {}
		self.id_to_coordinates = {}
		self.coordinates_to_id = {}
		self.link_capacity = 100 # in Gbps
		#print self.S
		# first generate all the switches
		self.create_switches()
		self.wire_network()
		return

	def create_switches(self):
		#num_switches = product = reduce((lambda x, y: x * y), self.S)
		#dims = [0] * self.L
		# First step is to create all of the switches, this can be done by just generating the cartesian product of all the coordinates
		all_coords = []
		for i in range(self.S[-1]):
			all_coords.append((i,))
		for dim2 in range(self.L - 2, -1, -1):
			coords1 = []
			for i in range(self.S[dim2]):
				coords1.append((i,))
			all_coords = cart_product(coords1, all_coords)
		swid = 0
		for coord in all_coords:
			self.adjacency_list[coord] = []
			self.id_to_coordinates[swid] = coord
			self.coordinates_to_id[coord] = swid
			swid += 1
		#print all_coords
		print "total number of switches = " + str(len(all_coords))
		return

	def wire_network(self):
		# Next step is to wire all of them together
		num_neighbors = 0
		for i in self.S:
			num_neighbors += (i - 1)
		for src in self.adjacency_list.keys():
			for dst in self.adjacency_list.keys():
				if (not self.share_dimension(src, dst)) or (dst in self.adjacency_list[src]):
					continue
				else:
					self.adjacency_list[src].append(dst)
			assert(len(self.adjacency_list[src]) == num_neighbors)
		return

	# checks to see if there is at least one dimension, returns True is so, and False otherwise
	def share_dimension(self, coord1, coord2):
		see_diff = False
		for i in range(len(coord1)):
			if coord1[i] == coord2[i]:
				continue
			elif coord1[i] != coord2[i] and see_diff:
				return False
			else:
				see_diff = True
		return see_diff


	def cheeger_constant(self):
		min_so_far = 10E10
		for i in range(1000):
			A, B = splitset(self.adjacency_list.keys())
			boundary_links = 0
			for elem in A:
				for neighbor in self.adjacency_list[elem]:
					if neighbor in B:
						boundary_links += 1
			min_so_far = min(float(boundary_links) / len(A), min_so_far)
		return min_so_far

	def adjacency_matrix(self):
		num_switches = len(self.adjacency_list.keys())
		mat = [0] * num_switches
		offset = 0
		coord_to_index = {}
		for sw_coord in self.adjacency_list.keys():
			coord_to_index[sw_coord] = offset
			mat[offset] = [0] * num_switches
			offset += 1
		for sw in self.adjacency_list.keys():
			for neighbor in self.adjacency_list[sw]:
				mat[coord_to_index[sw]][coord_to_index[neighbor]] += 1
		return mat

	def path_recur(self, remaining_entries):
		if len(remaining_entries) == 0:
			return []
		elif len(remaining_entries) == 1:
			return [remaining_entries]
		else:
			all_paths = []
			for i in range(len(remaining_entries)):
				copied_list = copy.copy(remaining_entries)
				copied_list.pop(i)
				collection = self.path_recur(copied_list)
				for item in collection:
					item.append(remaining_entries[i])
					all_paths.append(item)
			return all_paths

	# returns a set of paths (list of vertices which are id of switches) which are shortest paths 
	# between a src and dst
	def shortest_path_set(self, adj_matrix, src, dst):
		src_coord = self.id_to_coordinates[src]
		dst_coord = self.id_to_coordinates[dst]
		dims_diffs = []
		print src_coord
		print dst_coord
		for i in range(len(src_coord)):
			if src_coord[i] != dst_coord[i]:
				dims_diffs.append(i)
		all_dim_sequences = self.path_recur(dims_diffs)
		print all_dim_sequences
		all_paths = []
		# the path_recur function merely returns all permutations of sequences of dimensions to take
		# we still need to transform them into id's of the switches
		for dim_seq in all_dim_sequences:
			path = [src, ]
			curr_coord = list(self.id_to_coordinates[src])
			for dim_ind in dim_seq:				
				curr_coord[dim_ind] = dst_coord[dim_ind]				
				path.append(self.coordinates_to_id[tuple(curr_coord)])
			all_paths.append(path)
		return all_paths

	# routes a traffic matrix using wcmp
	def route_wcmp(traffic_matrix, network):
		adj_matrix = self.adjacency_matrix()
		num_switches = len(adj_matrix)
		traffic_load = [0] * num_switches # traffic load on each link between links
		for i in range(num_switches):
			traffic_load[i] = [0] * num_switches
		for src in range(num_switches):
			for dst in range(num_switches):
				if src == dst:
					continue
				path_total_weight = 0
				path_set = self.shortest_path_set(adj_matrix, src, dst)
				path_weight = []
				for path in path_set:	
					path_min_capacity = self.link_capacity
					curr_node = path[0]
					for i in range(1, len(path), 1):
						intermediate_swid = path[i]
						path_min_capacity = min(float(path_min_capacity), adj_matrix[curr_node][intermediate_swid] * self.link_capacity)
						curr_node = intermediate_swid
					path_weight.append(path_min_capacity)
					path_total_weight += path_min_capacity
				path_index = 0
				for path in path_set:
					curr_node = path[0]
					for i in range(1, len(path), 1):
						intermediate_swid = path[i]
						traffic_load[curr_node][intermediate_swid] += (traffic_matrix[src][dst] * (path_weight[path_index]/path_total_weight))
						curr_node = intermediate_swid
					path_index += 1
		# finally, compute the mlu
		mlu = 0
		for src in range(num_switches):
			for dst in range(num_switches):
				if src == dst:
					continue
				if adj_matrix[src][dst] <= 0:
					continue
				mlu = max(mlu, traffic_load[src][dst] / adj_matrix[src][dst] / self.link_capacity)
		return mlu
	# wire up network to optical switches within one dimension
	#def wire_network_intradimension():
	#	num_optical_switches = 0
	#	optical_switches = {}
	#	for dim in range(self.L):
	#		for ind in range(self.S[dim]):
				
	#	return

hx = HyperX(3, [3], 1, 1)
#hx = HyperX(4, [4], 1, 1)
adj_matrix = hx.adjacency_matrix()
path_set = hx.shortest_path_set(adj_matrix, 0, 26)
print path_set




	