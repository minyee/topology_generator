import sys, math
import hyperx_reconfigurable as rhpx
import hyperx as hpx
import traffic_generator as traffic_gen

taper = 0.5
dimensions = [3,3,3, 6]
tpx = rhpx.TaperedHyperX(len(dimensions), dimensions, 1, 1, taper)
#tpx.print_topology()
links_per_switch = int(taper * (dimensions[-1] - 1))
links_per_switch = max(links_per_switch, 1)
assert(tpx.check_topology_validity(links_per_switch))
# next generate the traffic matrices
tm = traffic_gen.generate_hpx_intergroup_bipartite(tpx)

mlu = tpx.route_wcmp(tm)
print "MLU with taper of : {} is: {}".format(taper, mlu)


taper = 1
tpx = rhpx.TaperedHyperX(len(dimensions), dimensions, 1, 1, taper)
#tpx.print_topology()
links_per_switch = int(taper * (dimensions[-1] - 1))
links_per_switch = max(links_per_switch, 1)
assert(tpx.check_topology_validity(links_per_switch))
# next generate the traffic matrices
tm = traffic_gen.generate_hpx_intergroup_bipartite(tpx)

mlu = tpx.route_wcmp(tm)
print "MLU with taper of : {} is: {}".format(taper, mlu)