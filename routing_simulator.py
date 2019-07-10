import sys, math
import hyperx_taperedstatic as thpx
import hyperx as hpx
import hyperx_reconfigurable as rhpx
import traffic_generator as traffic_gen
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import util

dimensions = [3,3,3, 6]
tpx = thpx.TaperedHyperX(len(dimensions), dimensions, 1, 1, 1.)
tm_cluster, interpod_tm_cluster = traffic_gen.generate_hpx_intergroup_clustered(tpx, 3)
total_switches = tpx.num_intragroup_switches() * tpx.num_groups()

trace_file_dir = "/Users/minyee/src/arpa_e/traces/"
trace_file = trace_file_dir + "milc_512"
real_traffic_pattern = traffic_gen.read_traffic_trace_files(trace_file)

real_traffic_pattern = traffic_gen.reshape_tm(real_traffic_pattern, total_switches)
## bipartite traffic
tm_bipartite, interpod_tm_bipartite = traffic_gen.generate_hpx_intergroup_bipartite(tpx)

assert(len(tm_bipartite) == len(real_traffic_pattern))
tm_bipartite = real_traffic_pattern
intergroup_tm = [0] * tpx.S[-1]
for group in range(tpx.S[-1]):
	intergroup_tm[group] = [0] * tpx.S[-1]
for i in range(len(tm_bipartite)):
	src_coord = tpx.id_to_coordinates[i]
	for j in range(len(tm_bipartite)):
		dst_coord = tpx.id_to_coordinates[j]
		intergroup_tm[src_coord[-1]][dst_coord[-1]] += float(tm_bipartite[i][j])
print intergroup_tm


tapering_factors = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.]
mlu_static = []
alu_static = []
mlu_reconf = []
alu_reconf = []
throughput_static = []
throughput_reconf = []
num_links_needed = []
#util.write_to_netbench_traffic_matrix(tm_bipartite, "/Users/minyee/src/netbench/temp/hyperX/topology_files/tm_bipartite.txt")
#util.write_to_netbench_traffic_matrix(tm_cluster, "/Users/minyee/src/netbench/temp/hyperX/topology_files/tm_cluster.txt")
#util.write_to_netbench_traffic_matrix(interpod_tm_bipartite, "/Users/minyee/src/netbench/temp/hyperX/topology_files/tm_interpod_bipartite.txt")
#util.write_to_netbench_traffic_matrix(interpod_tm_cluster, "/Users/minyee/src/netbench/temp/hyperX/topology_files/tm_interpod_cluster.txt")
system_size_static = []
system_size_reconf = []
for taper in tapering_factors:
	reconfigurable_hpx = rhpx.ReconfigurableHyperX(len(dimensions), dimensions, 1, 1, taper)
	reconfigurable_hpx.bandwidth_steer_ilp(taper, intergroup_tm)
	reconf_am = reconfigurable_hpx.get_adjacency_matrix()
	reconf_interpod_am = reconfigurable_hpx.get_interpod_adjacency_matrix()
	mlu_recon, alu_recon = reconfigurable_hpx.route_wcmp(tm_bipartite)
	taperedstatic_hpx = thpx.TaperedHyperX(len(dimensions), dimensions, 1, 1, taper)
	tapered_am = taperedstatic_hpx.get_adjacency_matrix()
	tapered_interpod_am = taperedstatic_hpx.get_interpod_adjacency_matrix()
	system_size_reconf.append(reconfigurable_hpx.get_num_servers(32))
	system_size_static.append(taperedstatic_hpx.get_num_servers(32))
	mlu_stat, alu_stat = taperedstatic_hpx.route_wcmp(tm_bipartite)
	mlu_static.append(mlu_stat)
	alu_static.append(alu_stat)
	mlu_reconf.append(mlu_recon)
	alu_reconf.append(alu_recon)
	throughput_reconf.append(1./mlu_recon)
	throughput_static.append(1./mlu_stat)
	num_links_needed.append(len(taperedstatic_hpx.adjacency_list[(0,0,0,0)]))
	tapered_string = "%1.1f" % taper
	tapered_string = tapered_string.replace(".", "p")
	#util.write_to_netbench_topology_format(tapered_am, "/Users/minyee/src/netbench/temp/hyperX/topology_files/tapered_" + tapered_string + ".topology")
	#util.write_to_netbench_topology_format(reconf_am, "/Users/minyee/src/netbench/temp/hyperX/topology_files/reconf_" + tapered_string + ".topology")
	#util.write_to_netbench_interpod_topology_format(tapered_interpod_am, "/Users/minyee/src/netbench/temp/hyperX/topology_files/tapered_interpod_" + tapered_string + ".topology")
	#util.write_to_netbench_interpod_topology_format(reconf_interpod_am, "/Users/minyee/src/netbench/temp/hyperX/topology_files/reconf_interpod_" + tapered_string + ".topology")

fig = plt.figure()
plt.plot(tapering_factors, system_size_static)
plt.plot(tapering_factors, system_size_static)
plt.ylabel('Num Servers')
plt.xlabel('Fraction of links left')
plt.legend(['Non-reconfigurable', 'Reconfigurable'])

fig = plt.figure()
plt.plot(tapering_factors, mlu_static)
plt.plot(tapering_factors, mlu_reconf)
plt.title('Maximum congestion')
plt.ylabel('Max congestion')
plt.xlabel('Fraction of links left')
plt.legend(['Non-reconfigurable', 'Reconfigurable'])

fig = plt.figure()
plt.plot(tapering_factors, alu_static)
plt.plot(tapering_factors, alu_reconf)
plt.title('Average congestion')
plt.ylabel('Ave congestion')
plt.xlabel('Fraction of links left')
plt.legend(['Non-reconfigurable', 'Reconfigurable'])

fig = plt.figure()
plt.plot(tapering_factors, num_links_needed)
plt.title('Number of ports per switch for network-facing connection')
plt.ylabel('Num Ports')
plt.xlabel('Fraction of links left')

fig = plt.figure()
plt.plot(tapering_factors, throughput_static)
plt.plot(tapering_factors, throughput_reconf)
plt.title('Throughput of Flows')
plt.ylabel('Throughput')
plt.xlabel('Fraction of links left')
plt.legend(['Non-Reconfigurable', 'Reconfigurable'])
plt.show()

## clustered traffic
intergroup_tm = [0] * tpx.S[-1]
for group in range(tpx.S[-1]):
	intergroup_tm[group] = [0] * tpx.S[-1]
for i in range(len(tm_cluster)):
	src_coord = tpx.id_to_coordinates[i]
	for j in range(len(tm_cluster)):
		dst_coord = tpx.id_to_coordinates[j]
		intergroup_tm[src_coord[-1]][dst_coord[-1]] += tm_cluster[i][j]


tapering_factors = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.]
mlu_static = []
alu_static = []
mlu_reconf = []
alu_reconf = []

for taper in tapering_factors:
	reconfigurable_hpx = rhpx.ReconfigurableHyperX(len(dimensions), dimensions, 1, 1, taper)
	reconfigurable_hpx.bandwidth_steer_ilp(taper, intergroup_tm)
	mlu_recon, alu_recon = reconfigurable_hpx.route_wcmp(tm_cluster)
	taperedstatic_hpx = thpx.TaperedHyperX(len(dimensions), dimensions, 1, 1, taper)
	mlu_stat, alu_stat = taperedstatic_hpx.route_wcmp(tm_cluster)
	mlu_static.append(mlu_stat)
	alu_static.append(alu_stat)
	mlu_reconf.append(mlu_recon)
	alu_reconf.append(alu_recon)


fig = plt.figure()
plt.plot(tapering_factors, mlu_static)
plt.plot(tapering_factors, mlu_reconf)
plt.title('Maximum congestion')
plt.ylabel('Max congestion')
plt.xlabel('Fraction of links left')
plt.legend(['Non-reconfigurable', 'Reconfigurable'])

fig = plt.figure()
plt.plot(tapering_factors, alu_static)
plt.plot(tapering_factors, alu_reconf)
plt.title('Average congestion')
plt.ylabel('Ave congestion')
plt.xlabel('Fraction of links left')
plt.legend(['Non-reconfigurable', 'Reconfigurable'])

plt.show()	