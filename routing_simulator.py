import sys, math
import hyperx_taperedstatic as thpx
import hyperx as hpx
import hyperx_reconfigurable as rhpx
import traffic_generator as traffic_gen

taper = 0.5
dimensions = [3,3,3, 6]
tpx = thpx.TaperedHyperX(len(dimensions), dimensions, 1, 1, taper)
#tpx.print_topology()
links_per_switch = int(taper * (dimensions[-1] - 1))
links_per_switch = max(links_per_switch, 1)
assert(tpx.check_topology_validity(links_per_switch))
# next generate the traffic matrices
tm = traffic_gen.generate_hpx_intergroup_bipartite(tpx)

mlu = tpx.route_wcmp(tm)
print "MLU with taper of : {} is: {}".format(taper, mlu)

taper = 0.5 
tpx = rhpx.ReconfigurableHyperX(len(dimensions), dimensions, 1, 1, taper)
intergroup_tm = [0] * tpx.S[-1]
for group in range(tpx.S[-1]):
	intergroup_tm[group] = [0] * tpx.S[-1]
for i in range(len(tm)):
	src_coord = tpx.id_to_coordinates[i]
	for j in range(len(tm)):
		dst_coord = tpx.id_to_coordinates[j]
		intergroup_tm[src_coord[-1]][dst_coord[-1]] += tm[i][j]
tpx.bandwidth_steer_ilp(taper, intergroup_tm)
mlu = tpx.route_wcmp(tm)
print "MLU with taper of : {} and bandwidth steered is: {}".format(taper, mlu)

taper = 1
tpx = thpx.TaperedHyperX(len(dimensions), dimensions, 1, 1, taper)
#tpx.print_topology()
links_per_switch = int(taper * (dimensions[-1] - 1))
links_per_switch = max(links_per_switch, 1)
assert(tpx.check_topology_validity(links_per_switch))
# next generate the traffic matrices
tm = traffic_gen.generate_hpx_intergroup_bipartite(tpx)

mlu = tpx.route_wcmp(tm)
print "MLU with taper of : {} is: {}".format(taper, mlu)