import sys
from functools import reduce
import random

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
		print self.S
		# first generate all the switches
		self.create_switches()

		self.wire_network()
		return

	# Given a switch id, returns the high-dimensional coordinates of the switch
	def id_to_coordinates(self, id):
		return self.id_to_indices[id]

	def coordinates_to_id(self, indices):
		offset = 0
		i = 0
		for dim in range(self.L - 1, 0, -1):
			mult = 1
			for dim2 in range(dim - 1, 0, -1):
				mult *= self.S[dim2]
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
		for coord in all_coords:
			self.adjacency_list[coord] = []
		print all_coords
		print "total number of switches = " + str(len(all_coords))
		return

	def wire_network(self):
		# Next step is to wire all of them together
		num_neighbors = 1
		for i in self.S:
			num_neighbors *= (i - 1)
		for src in self.adjacency_list.keys():
			for dst in self.adjacency_list.keys():
				if src == dst:
					continue
				if self.share_dimension(src, dst) or (dst in self.adjacency_list[src]):
					continue
				self.adjacency_list[src].append(dst)
			assert(len(self.adjacency_list[src]) == num_neighbors)
		return

	# checks to see if there is at least one dimension, returns True is so, and False otherwise
	def share_dimension(self, coord1, coord2):
		for i in range(len(coord1)):
			if coord1[i] == coord2[i]:
				return True
		return False


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


	# wire up network to optical switches within one dimension
	#def wire_network_intradimension():
	#	num_optical_switches = 0
	#	optical_switches = {}
	#	for dim in range(self.L):
	#		for ind in range(self.S[dim]):
				
	#	return

#hx = HyperX(3, [3], 1, 1)
#hx = HyperX(4, [4], 1, 1)
#print "cheeger constant -> {}".format(hx.cheeger_constant())




	