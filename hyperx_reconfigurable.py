import sys
from gurobipy import *
import random
import copy
import math
import networkx as nx
import util 
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

class ReconfigurableHyperX:
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

	def compute_logical_intergroup_connectivity(self, intergroup_traffic_matrix, intergroup_links):
		logical_intergroup_topology = [0] * self.S[-1]
		for group in range(self.S[-1]):
			logical_intergroup_topology[group] = [0] * self.S[-1]
		# first we need to figure out what the optimal logical intergroup topology is 
		model = Model("Logical topology lp")
		decision_vars = {}
		intergroup_traffic_matrix = util.normalize_matrix(intergroup_traffic_matrix)
		for i in range(self.S[-1] - 1):
			for j in range(i, self.S[-1], 1):
				decision_vars[(i, j)] = model.addVar(0,intergroup_links, 0., GRB.CONTINUOUS, "{}".format((i,j)))
		mu = model.addVar(0, GRB.INFINITY, 1., GRB.CONTINUOUS, "mu")
		# ensuring the mu is the maximum traffic scale up1
		for i in range(self.S[-1] - 1):
			for j in range(i + 1, self.S[-1], 1):
				model.addConstr(lhs=decision_vars[(i,j)], sense=GRB.GREATER_EQUAL, rhs=mu * (intergroup_traffic_matrix[i][j] + intergroup_traffic_matrix[j][i]))
		# number of incoming and outgoing links constraint
		for curr_group in range(self.S[-1]):
			link_constraint = LinExpr()
			for target_group in range(self.S[-1]):
				if curr_group != target_group:
					ind = (min(curr_group, target_group), max(curr_group, target_group))
					link_constraint.addTerms(1., decision_vars[(ind)])
			model.addConstr(lhs=link_constraint, sense=GRB.EQUAL, rhs=intergroup_links)
		model.setObjective(mu, GRB.MAXIMIZE)
		model.optimize()
		for i in range(self.S[-1] - 1):
			for j in range(i, self.S[-1], 1):
				var = model.getVarByName("{}".format((i,j)))
				logical_intergroup_topology[i][j] = var.x
				logical_intergroup_topology[j][i] = var.x # check if this should be populated here
		print "\n\n\n\nintergroup_traffic_matrix\n\n\n\n"
		self.print_matrix(intergroup_traffic_matrix)
		print "\n\n\n\nlogical_intergroup_topology\n\n\n\n"
		self.print_matrix(logical_intergroup_topology)
		return logical_intergroup_topology

	def bandwidth_steer_ilp(self, taper, intergroup_traffic_matrix):
		links_per_switch = int(taper * (self.S[-1] - 1))
		links_per_switch = max(links_per_switch, 1)
		decision_vars = {}
		num_intra_group_switches = 1
		num_intragroup_links = 0
		for dim in range(len(self.S) - 1):
			num_intra_group_switches *= self.S[dim]
			num_intragroup_links += (self.S[dim] - 1) 		
		num_intergroup_links = num_intra_group_switches * links_per_switch
		logical_intergroup_topology = self.compute_logical_intergroup_connectivity(intergroup_traffic_matrix, num_intergroup_links)
		model = Model("topology ilp")
		# form all the switch coordinates
		switch_in_group = [0] * self.S[-1]
		for group in range(self.S[-1]):
			coord = [0] * (len(self.S) - 1)
			switch_in_group[group] = []
			for intra_group_switch in range(num_intra_group_switches):
				switch_in_group[group].append(list(coord))
				coord = self.increment_switch_coord(coord, self.S[:-1])
		# prelude, set up all the decision optimization variables first
		connectivity_variables = {}
		for src_group in range(self.S[-1] - 1):
			for dst_group in range(src_group + 1, self.S[-1], 1):
				for coord1 in switch_in_group[src_group]:
					src_coord = tuple(coord1 + [src_group,])
					for coord2 in switch_in_group[dst_group]:
						dst_coord = tuple(coord2 + [dst_group,])
						key = (self.coordinates_to_id[src_coord], self.coordinates_to_id[dst_coord])
						connectivity_variables[key] = model.addVar(0, 1, 0., GRB.INTEGER, "(" + str(self.coordinates_to_id[src_coord]) + ", " + str(self.coordinates_to_id[dst_coord]) + ")")

		# now add in the the connectivity constraints
		intergroup_connectivity_constraints = {}
		for src_group in range(self.S[-1] - 1):
			for dst_group in range(src_group + 1, self.S[-1], 1):
				intergroup_connectivity_constraints[(src_group, dst_group)] = LinExpr()				
				for coord1 in switch_in_group[src_group]:
					src_coord = tuple(coord1 + [src_group])
					for coord2 in switch_in_group[dst_group]:
						dst_coord = tuple(coord2 + [dst_group])
						key = (self.coordinates_to_id[src_coord], self.coordinates_to_id[dst_coord])
						intergroup_connectivity_constraints[(src_group, dst_group)].addTerms(1., connectivity_variables[key])
		# figure out where to add the constraints here
		for src_group in range(self.S[-1] - 1):
			for dst_group in range(src_group + 1, self.S[-1], 1):
				model.addConstr(lhs=intergroup_connectivity_constraints[(src_group, dst_group)], sense=GRB.GREATER_EQUAL, rhs=math.floor(logical_intergroup_topology[src_group][dst_group]), name="connectivity_lb, {} - {}".format(src_group, dst_group))
				model.addConstr(lhs=intergroup_connectivity_constraints[(src_group, dst_group)], sense=GRB.LESS_EQUAL, rhs=math.ceil(logical_intergroup_topology[src_group][dst_group]), name="connectivity_ub, {} - {}".format(src_group, dst_group))
		# Finally, optimize the model, and if successful, then create the topology
		try:
			topology_matrix = [0] * self.S[-1]
			for i in range(self.S[-1]):
				topology_matrix[i] = [0] * self.S[-1]
			model.optimize()
			for src_group in range(self.S[-1] - 1):
				for dst_group in range(src_group + 1, self.S[-1], 1):
					for coord1 in switch_in_group[src_group]:
						src_coord = tuple(coord1 + [src_group])
						for coord2 in switch_in_group[dst_group]:
							dst_coord = tuple(coord2 + [dst_group])
							var_name = "(" + str(self.coordinates_to_id[src_coord]) + ", " + str(self.coordinates_to_id[dst_coord]) + ")"
							var = model.getVarByName(var_name)
							if var.x >= 1:
								topology_matrix[src_group][dst_group] += var.x
								topology_matrix[dst_group][src_group] += var.x
								self.adjacency_list[src_coord].append(dst_coord)
								self.adjacency_list[dst_coord].append(src_coord)
			for src_group in range(self.S[-1]):
				row = "| "
				for entry in topology_matrix[src_group]:
					row += (str(entry) + " ")
				row += "|\n"
				print row
		except GurobiError as e:
			print ("Error code " + str(e. errno ) + ": " + str(e))
		except AttributeError :
			print ("Encountered an attribute error ")
		return 

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
		for src_coord in self.adjacency_list.keys():
			src = self.coordinates_to_id[src_coord]
			mat[src] = [0] * num_switches
			for dst_coord in self.adjacency_list[src_coord]:
				dst = self.coordinates_to_id[dst_coord]
				mat[src][dst] += 1
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
	def shortest_path_set(self, adj_matrix, src, dst, G):
		# sense if the src and dst are separated in the tapered dimension. That is, if the last coordinate is the same
		# if they are the same, can do the normal way of using shortest_path_set
		# if not, then need to find an entrance switch that has a direct connectivity to the final dimension
		src_coord = self.id_to_coordinates[src]
		dst_coord = self.id_to_coordinates[dst]
		dims_diffs = []
		for i in range(len(src_coord)):
			if src_coord[i] != dst_coord[i]:
				dims_diffs.append(i)
		if self.L - 1 not in dims_diffs:
			all_dim_sequences = self.path_recur(dims_diffs)
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
		else:
			# first, form the graph object
			# now that the graph has been built fully, time to return all paths
			all_paths = []
			for path in nx.all_shortest_paths(G, source=src, target=dst):
				all_paths.append(path)
			return all_paths


	# routes a traffic matrix using wcmp
	def route_wcmp(self, traffic_matrix):
		adj_matrix = self.adjacency_matrix()
		num_switches = len(adj_matrix)
		traffic_load = [0] * num_switches # traffic load on each link between links
		G = nx.Graph()
		G.add_nodes_from(range(self.num_groups() * self.num_intragroup_switches()))
		for i in range(len(adj_matrix)):
			for j in range(len(adj_matrix)):
				if i != j and adj_matrix[i][j] > 0:
					G.add_edge(i, j)
		for i in range(num_switches):
			traffic_load[i] = [0] * num_switches
		for src in range(num_switches):
			for dst in range(num_switches):
				if src != dst:
					path_total_weight = 0
					path_set = self.shortest_path_set(adj_matrix, src, dst, G)
					path_weight = []
					#print "src: {} coord - {}, dst: {} coord - {}".format(src, self.id_to_coordinates[src], dst, self.id_to_coordinates[dst])
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
		mlu = 0.
		alu = 0.
		link_count = 0
		for src in range(num_switches):
			for dst in range(num_switches):
				if src == dst or adj_matrix[src][dst] <= 0:
					continue
				else:
					mlu = max(mlu, traffic_load[src][dst] / adj_matrix[src][dst] / self.link_capacity)
					alu += (traffic_load[src][dst] / adj_matrix[src][dst] / self.link_capacity)
					link_count += adj_matrix[src][dst]
		return mlu, alu / link_count




	def print_topology(self):
		for coord in self.adjacency_list.keys():
			print "\n\nCurrent coord: {}, id: {}".format(coord, self.coordinates_to_id[coord])
			for neighbors in self.adjacency_list[coord]:
				print "coord: {}, id, {}".format(neighbors, self.coordinates_to_id[neighbors])

	def print_matrix(self, matrix):
		for row in matrix:
			string = "| "
			for elem in row:
				string += (str(elem) + ", ")
			string += "|"
			print string


	def num_groups(self):
		return self.S[-1]

	# number of switches in all bar the final dimensions
	def num_intragroup_switches(self):
		num = 1
		for i in self.S[:-1]:
			num *= i
		return num

	def get_adjacency_matrix(self):
		num_switches = len(self.adjacency_list.keys())
		adj_matrix = [0] * num_switches
		for src in self.adjacency_list.keys():
			src_id = self.coordinates_to_id[src]
			adj_matrix[src_id] = [0] * num_switches
			for neighbor in self.adjacency_list[src]:
				neighbor_id = self.coordinates_to_id[neighbor]
				adj_matrix[src_id][neighbor_id] += 1
		return adj_matrix

	def get_interpod_adjacency_matrix(self):
		interpod_matrix = [0] * self.S[-1]
		for i in range(self.S[-1]):
			interpod_matrix[i] = [0] * self.S[-1]
		for src in self.adjacency_list.keys():
			src_pod = src[-1]
			for dst in self.adjacency_list[src]:
				dst_pod = dst[-1]
				if src_pod != dst_pod:
					interpod_matrix[src_pod][dst_pod] += 1
		return interpod_matrix

	def get_num_servers(self, eps_radix):
		num_servers = 0
		for switch in self.adjacency_list.keys():
			network_ports = len(self.adjacency_list[switch])
			num_servers += (eps_radix - network_ports)
		return num_servers
#thx = TaperedHyperX(4, [3,3,3,6], 1, 1, 0.5)
#adj_mat = thx.adjacency_matrix()


	