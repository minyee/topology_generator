import sys
import copy

def normalize_matrix(matrix):
	agg = 0.
	entries = len(matrix)
	for i in range(entries):
		for j in range(entries):
			agg += matrix[i][j]
	matrix_copy = copy.copy(matrix)
	for i in range(entries):
		for j in range(entries):
			matrix_copy[i][j] = matrix[i][j] / agg
	return matrix_copy


def write_to_netbench_topology_format(adj_matrix, output_filename):
	num_switches = len(adj_matrix)
	num_edges = 0
	for i in range(num_switches):
		for j in range(num_switches):
			num_edges += adj_matrix[i][j]
	str_builder = "\n"
	str_builder += ("|V|=" + str(num_switches) + "\n")
	str_builder += ("|E|=" + str(num_edges) + "\n")
	str_builder += ("ToRs=incl_range(0," + str(num_switches - 1) + ")\n")
	str_builder += ("Servers=incl_range(0," + str(num_switches - 1) + ")\n")
	str_builder += ("Switches=set()\n\n# Links\n")
	for i in range(num_switches):
		for j in range(num_switches):
			num_links = adj_matrix[i][j]
			while (num_links > 0):
				str_builder += (str(i) + " " + str(j) + "\n")
				num_links -= 1
	with open(output_filename, "w+") as f:
		f.write(str_builder)
	return

def write_to_netbench_interpod_topology_format(interpod_matrix, output_filename):
	num_switches = len(interpod_matrix)
	num_edges = 0
	for i in range(num_switches):
		for j in range(num_switches):
			num_edges += interpod_matrix[i][j]
	str_builder = "\n"
	str_builder += ("|V|=" + str(num_switches) + "\n")
	str_builder += ("|E|=" + str(num_edges) + "\n")
	str_builder += ("ToRs=incl_range(0," + str(num_switches - 1) + ")\n")
	str_builder += ("Servers=incl_range(0," + str(num_switches - 1) + ")\n")
	str_builder += ("Switches=set()\n\n# Links\n")
	for i in range(num_switches):
		for j in range(num_switches):
			num_links = interpod_matrix[i][j]
			while (num_links > 0):
				str_builder += (str(i) + " " + str(j) + "\n")
				num_links -= 1
	with open(output_filename, "w+") as f:
		f.write(str_builder)
	return

def write_to_netbench_traffic_matrix(traffic_matrix, output_filename):
	str_builder = "#tor_pair_id,src,dst,pdf_num_bytes\n"
	pair_id = 0
	num_switches = len(traffic_matrix)
	norm_tm = normalize_matrix(traffic_matrix)
	num_switches = len(traffic_matrix)
	for i in range(num_switches):
		for j in range(num_switches):
			if i != j and norm_tm[i][j] > 0.:
				str_builder += ("" + str(pair_id) + "," + str(num_switches + i) + "," + str(num_switches + j) + "," + ("%.6g" % norm_tm[i][j]) + "\n")
				pair_id += 1
	str_builder += "\n"
	with open(output_filename, "w+") as f:
		f.write(str_builder)
	return
