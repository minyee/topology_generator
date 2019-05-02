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

class TaperedHyperX:
	def __init__(self, L_arg, S_arg, K_arg, T_arg, taper):
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
		self.wire_static_network()
		self.wire_tapered_dimension(taper)
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

	# in this case, we first don't connect the final dimension
	def wire_static_network(self):
		num_neighbors = 0
		for i in self.S[:-1]:
			num_neighbors += (i - 1)
		for src in self.adjacency_list.keys():
			for dst in self.adjacency_list.keys():
				if (not self.share_dimension(src, dst)) or (dst in self.adjacency_list[src]):
					continue
				else:
					self.adjacency_list[src].append(dst)
			assert(len(self.adjacency_list[src]) == num_neighbors)
		return

	# increments the coordinates of the switch, dimension_size is the number of indices max in each dimension
	def increment_switch_coord(self, coord, dimension_size):
		assert(len(dimension_size) == len(coord))
		coord[-1] = (coord[-1] + 1) % dimension_size[-1]
		if coord[-1] != 0:
			return coord
		else:
			carry_forward = True
			for i in range(len(dimension_size) - 2, -1, -1):				
				if not carry_forward:
					break
				else:
					coord[i] = (coord[i] + 1) % self.S[i]
					if coord[i] != 0:
						carry_forward = False
			return coord

	def find_minimally_connected_group(self, intergroup_connectivity, curr_group):
		minimum = sys.maxsize
		min_group = -1
		for group in range(len(intergroup_connectivity)):
			if curr_group == group:
				continue
			else:
				if minimum > intergroup_connectivity[group]:
					minimum = intergroup_connectivity[group]
					min_group = group
		return min_group

	# selects the current group to connect links, returns true if we can continue wiring,
	# returns false otherwise
	def select_source_group(self, formed_links, maxlink):
		idx = -1
		minimum = sys.maxsize
		min_idx = 0
		continue_wiring = False
		wirable_groups = 0
		print formed_links
		for entry in formed_links:
			idx += 1
			if entry < maxlink:
				wirable_groups += 1
			if minimum > entry:
				minimum = entry
				min_idx = idx
		print min_idx
		return (wirable_groups >= 2), min_idx


	#def wire_tapered_dimension(self, taper):
	#	links_per_switch = int(taper * (self.S[-1] - 1))
	#	links_per_switch = max(links_per_switch, 1) # you need at least two links to ensure full connectivity, in which case the final dimension is just a ring
	#	group_to_group_links_formed = [0] * self.S[-1] # number of connections formed between a group to another group
	#	switch_in_group = [0] * self.S[-1] # coordinates of switches in each group
	#	num_intra_group_switches = 1
	#	num_intragroup_links = 0
	#	for dim in range(len(self.S) - 1):
	#		num_intra_group_switches *= self.S[dim]
	#		num_intragroup_links += (self.S[dim] - 1)
	#	if links_per_switch * num_intra_group_switches % 2 == 1:
	#		links_per_switch += 1
	#		links_per_switch = min(links_per_switch, self.S[-1] - 1)
	#	for i in range(self.S[-1]):
	#		group_to_group_links_formed[i] = [0] * self.S[-1]
	#		switch_in_group[i] = []
	#		coord = [0] * (len(self.S) - 1)
	#		for intra_group_switch in range(num_intra_group_switches):
	#			switch_in_group[i].append(list(coord))
	#			coord = self.increment_switch_coord(coord, self.S[:-1])
	#	wires_used = [0] * self.S[-1]
	#	max_link_per_group = links_per_switch * num_intra_group_switches
	#	continue_wiring = True
	#	src_group = 0
	#	while continue_wiring:	
	#		print "src_group: {}".format(src_group)
	#		src_switch = switch_in_group[src_group][-1] + [src_group,]
	#		switch_in_group[src_group].pop()
	#		while (num_intragroup_links + links_per_switch) == len(self.adjacency_list[tuple(src_switch)]):
	#				src_switch =  switch_in_group[target_group].pop()
	#				switch_in_group[src_group].pop()
	#		for link in range((num_intragroup_links + links_per_switch) - len(self.adjacency_list[tuple(src_switch)])):	
	#			continue_finding, target_group = self.find_minimally_connected_group(group_to_group_links_formed[src_group], src_group)
	#			while not continue_finding:
	#
	#				continue_finding, target_group = self.find_minimally_connected_group(group_to_group_links_formed[src_group], src_group)
	#			if len(switch_in_group[target_group]) == 0:
	#				print "intergroup_connectivity \n {}".format(group_to_group_links_formed) 
	#			target_switch = tuple(switch_in_group[target_group][-1]) + (target_group,)
	#			group_to_group_links_formed[src_group][target_group] += 1
	#			group_to_group_links_formed[target_group][src_group] += 1
	#			self.adjacency_list[tuple(src_switch)].append(target_switch)
	#			self.adjacency_list[target_switch].append(tuple(src_switch))
	#			print "wtf\n\n"
	#			wires_used[src_group] += 1
	#			wires_used[target_group] += 1
	#			if num_intragroup_links + links_per_switch == len(self.adjacency_list[target_switch]):
	#				switch_in_group[target_group].pop()
	#		print wires_used
	#		continue_wiring, src_group = self.select_source_group(wires_used, max_link_per_group)
	#	return

	def wire_tapered_dimension(self, taper):
		links_per_switch = int(taper * (self.S[-1] - 1))
		# you need at least two links to ensure full connectivity, in which case the final dimension is just a ring
		links_per_switch = max(links_per_switch, 1)
		num_intra_group_switches = 1
		num_intragroup_links = 0
		for dim in range(len(self.S) - 1):
			num_intra_group_switches *= self.S[dim]
			num_intragroup_links += (self.S[dim] - 1) 		

		# Step 1: figure out interblock link count first
		# the number of links connectted in a full mesh
		l_ij = int(num_intra_group_switches * links_per_switch / (self.S[-1] - 1))
		intergroup_topology = [0] * self.S[-1]
		links_left = [0] * self.S[-1]
		switch_in_group = [0] * self.S[-1]
		for i in range(self.S[-1]):
			intergroup_topology[i] = [l_ij] * self.S[-1]
			intergroup_topology[i][i] = 0
			links_left[i] = links_per_switch * num_intra_group_switches - (l_ij * (self.S[-1] - 1))
			switch_in_group[i] = []
			coord = [0] * (len(self.S) - 1)
			for intra_group_switch in range(num_intra_group_switches):
				switch_in_group[i].append(list(coord))
				coord = self.increment_switch_coord(coord, self.S[:-1])
		# now figure out the remainder
		for i in range(self.S[-1]):
			for j in range(self.S[-1]):
				if links_left[i] == 0:
					break
				if links_left[j] == 0:
					continue
				if i == j or intergroup_topology[i][j] > l_ij:
					continue
				else:
					intergroup_topology[i][j] += 1
					intergroup_topology[j][i] += 1
					links_left[i] -= 1
					links_left[j] -= 1
		total_links = num_intragroup_links + links_per_switch
		print "\n\ntotal links: {}\n\n".format(total_links)
		# Step 2: next figure out to distribute each blocks links amongst the switches
		for coord in switch_in_group[0]:
			for src_group in range(self.S[-1]):
				src_coord = tuple(coord + [src_group,])
				src_links_left = total_links - len(self.adjacency_list[src_coord])
				#print "\n\n entering coord: {}".format(src_coord)
				for dst_group in range(src_group + 1, self.S[-1], 1):
					print (src_group, dst_group)
					dst_coord = tuple(coord + [dst_group,])
					dst_links_left = total_links - len(self.adjacency_list[dst_coord])
					if intergroup_topology[src_group][dst_group] == 0 or dst_links_left == 0:
						continue
					elif src_links_left == 0:
						break
					else:
						#print "dst_group: {} ".format(dst_group)
						self.adjacency_list[src_coord].append(dst_coord)
						self.adjacency_list[dst_coord].append(src_coord)
						intergroup_topology[src_group][dst_group] -= 1
						intergroup_topology[dst_group][src_group] -= 1
						src_links_left -= 1
						dst_links_left -= 1
				#if src_links_left > 0:
					#print "coord is: {} and there are {} links left\n\n".format(coord + [src_group], src_links_left)
				#assert(src_links_left == 0)
		#print intergroup_topology
		return
	def wire_tapered_dimension_ilp(self, taper):
		links_per_switch = int(taper * (self.S[-1] - 1))
		links_per_switch = max(links_per_switch, 1)
		decision_vars = {}
		num_intra_group_switches = 1
		num_intragroup_links = 0
		for dim in range(len(self.S) - 1):
			num_intra_group_switches *= self.S[dim]
			num_intragroup_links += (self.S[dim] - 1) 		
		# add the radix constraints of each group
		for group in range(self.S[-1]):
			decision_vars[group] = {}
			coord = [0] * (len(self.S) - 1)
			# constraint_obj = 
			for intra_group_switch in range(num_intra_group_switches):
				switch_in_group[i].append(list(coord))
				coord = tuple(self.increment_switch_coord(coord, self.S[:-1]) + [group,])
				decision_vars[group][coord] = model.addVar(0, links_per_switch, 0., GRB.INTEGER, str(coord))
				+= decision_vars

		# add the radix constraints

		# add the connectivity constraints
		# finally add the symmetric connectivity matrix constraint

	# check validity of topology
	def check_topology_validity(self, final_dim_link):
		supposed_num_neighbors = 0
		for i in self.S:
			supposed_num_neighbors += (i - 1)
		supposed_num_neighbors = supposed_num_neighbors - (self.S[-1] - 1) + final_dim_link
		for coord in self.adjacency_list.keys():
			seen = {}
			num_neighbors = 0
			for neighbor in self.adjacency_list[coord]:
				if neighbor in seen.keys():
					return False
				else:
					seen[neighbor] = True
					num_neighbors += 1
			#assert(num_neighbors == supposed_num_neighbors)
		return True

	#def wire_tapered_dimension(self, taper):
	#	links_per_switch = int(taper * (self.S[-1] - 1))
	#	links_per_switch = max(links_per_switch, 1) # you need at least two links to ensure full connectivity, in which case the final dimension is just a ring
	#	group_to_group_links_formed = [0] * self.S[-1] # number of connections formed between a group to another group
	#	switch_in_group = [0] * self.S[-1] # coordinates of switches in each group
	#	num_intra_group_switches = 1
	#	num_intragroup_links = 0
	#	print self.S
	#	for dim in range(len(self.S) - 1):
	#		num_intra_group_switches *= self.S[dim]
	#		num_intragroup_links += (self.S[dim] - 1)
	#	if links_per_switch * num_intra_group_switches % 2 == 1:
	#		links_per_switch += 1
	#		links_per_switch = min(links_per_switch, self.S[-1] - 1)
	#	for i in range(self.S[-1]):
	#		group_to_group_links_formed[i] = [0] * self.S[-1]
	#		switch_in_group[i] = []
	#		coord = [0] * (len(self.S) - 1)
	#		for intra_group_switch in range(num_intra_group_switches):
	#			switch_in_group[i].append(list(coord))
	#			coord = self.increment_switch_coord(coord, self.S[:-1])
	#	print group_to_group_links_formed
	#	for src_group in range(self.S[-1]):
	#		while len(switch_in_group[src_group]) > 0:
	#			src_switch = switch_in_group[src_group].pop() + [src_group,]								
	#			for link in range((num_intragroup_links + links_per_switch) - len(self.adjacency_list[tuple(src_switch)])):	
	#				target_group = self.find_minimally_connected_group(group_to_group_links_formed[src_group], src_group)
	#				# first, try this
	#				if len(switch_in_group[target_group]) == 0:
	#					print "intergroup_connectivity \n {}".format(group_to_group_links_formed) 
	#				target_switch = tuple(switch_in_group[target_group][-1]) + (target_group,)
	#				group_to_group_links_formed[src_group][target_group] += 1
	#				group_to_group_links_formed[target_group][src_group] += 1
	#				self.adjacency_list[tuple(src_switch)].append(target_switch)
	#				self.adjacency_list[target_switch].append(tuple(src_switch))
	#				if num_intragroup_links + links_per_switch == len(self.adjacency_list[target_switch]):
	#					switch_in_group[target_group].pop()				
	#	return


	# checks to see if there is at least one dimension, returns True is so, and False otherwise
	def share_dimension(self, coord1, coord2):
		diff_sum = 0
		# connect everything but the final dimension
		for i in range(len(coord1)):
			if coord1[i] != coord2[i]:
				if diff_sum == 0 and (i == (len(coord1) - 1)):
					continue
				else:
					diff_sum += 1
		return (diff_sum == 1)



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

	def print_topology(self):
		for coord in self.adjacency_list.keys():
			print "\n\nCurrent coord: {}, id: {}".format(coord, self.coordinates_to_id[coord])
			for neighbors in self.adjacency_list[coord]:
				print "coord: {}, id, {}".format(neighbors, self.coordinates_to_id[neighbors])

	def num_groups(self):
		return self.S[-1]

	# number of switches in all bar the final dimensions
	def num_intragroup_switches(self):
		num = 1
		for i in self.S[:-1]:
			num *= i
		return num

#thx = TaperedHyperX(4, [3,3,3,6], 1, 1, 0.5)
#adj_mat = thx.adjacency_matrix()
#thx.print_topology()


	