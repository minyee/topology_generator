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
	def coonect_src_dst_switches(self, src, dst):
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
				self.coonect_src_dst_switches(src, (group_id * self.a) + offset)
			# next, connect to one other switch
			curr_offset = src % self.a
			target_group = (group_id + (curr_offset + 1)) % self.g
			self.coonect_src_dst_switches(src, target_group * self.a + curr_offset)
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

#dfly = CanonicalDragonfly(5, 2)