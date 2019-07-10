import numpy as np
import math

def is_d_regular(adjacency_matrix):
	found_d = False
	d_reg = 0
	for i in range(len(adjacency_matrix)):
		num_edge = 0
		for neighbor in adjacency_matrix[i]:
			num_edge += neighbor
		if not found_d:
			found_d = True
			d_reg = num_edge
		elif num_edge != d_reg:
			return False, -1
	return True, d_reg




# cheeger bound for d-regular graphs
def cheeger_bound(adjacency_matrix):
	boolean, d = is_d_regular(adjacency_matrix)

	if not boolean:
		return (-1, -1)
	adjacency = np.array(adjacency_matrix)
	
	w, v = np.linalg.eig(adjacency)
	w = sorted(w, reverse=True)
	lambda1 = w[1]
	lb = 0.5 * (d - lambda1)
	ub = math.sqrt(2 * d * (d - lambda1))
	return (lb, ub)

