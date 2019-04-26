import sys
from functools import reduce
import random
import hyperx


# For now can only support canonical dragonfly, where a = g - 1, h = 1, 
class CanonicalDragonfly:
	def __init__(self, g, p):
		self.g = g
		self.a = g - 1
		self.p = p # concentration factor of the dragonfly topology
		self.adjacency_list = {}
		self.build_all_switches()
		self.connect_all_switches()

	def build_all_switches(self):
		index = 0
		num_switches = self.g * self.a
		for i in range(num_switches):
			self.adjacency_list[i] = []

	# given an id of a switch, returns the group number of the switch
	def group_id(self, swid):
		return swid / self.a

	def intra_group_id():
		return swid % self.a

	# forms a directed link from src switch to dest switch
	# no error checks
	def connect_src_dst_switches(self, src, dst):
		self.adjacency_list[src].append(dst)

	def connect_all_switches(self):
		for src in self.adjacency_list.keys():
			group_id = self.group_id(src)
			# first form the internal group connections
			intra_group_neighbors = range(self.a)
			for offset in intra_group_neighbors:
				if (group_id * self.a) + offset == src:
					continue
				# otherwise, form a directed link\
				self.connect_src_dst_switches(src, (group_id * self.a) + offset)
		group_offset = [0] * self.g
		for src_group in range(self.g - 1):
			for target_group in range(src_group + 1, self.g, 1):
				sw1 = src_group * self.a + group_offset[src_group]
				sw2 = target_group * self.a + group_offset[target_group]
				self.connect_src_dst_switches(sw1, sw2)
				self.connect_src_dst_switches(sw2, sw1)
				group_offset[src_group] += 1
				group_offset[target_group] += 1
		return

	# computes the expandability of the network graph
	def cheeger_constant(self):
		trials = 2000
		cc = [0] * trials
		for i in range(trials):
			A, B = hyperx.splitset(self.adjacency_list.keys())
			boundary_links = 0
			for elem in A:
				for neighbor in self.adjacency_list[elem]:
					if neighbor in B:
						boundary_links += 1
			cc[i] = (float(boundary_links) / len(A))
		return min(cc)

	def adjacency_matrix(self):
		num_switches = len(self.adjacency_list.keys())
		mat = [0] * num_switches
		for i in range(num_switches):
			mat[i] = [0] * num_switches
		for sw in self.adjacency_list.keys():
			for neighbor in self.adjacency_list[sw]:
				mat[sw][neighbor] += 1
		return mat

#dfly = CanonicalDragonfly(5, 2)