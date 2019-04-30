import sys, math
import hyperx_reconfigurable as rhpx
import hyperx as hpx
import traffic_generator as traffic_gen

taper = 0.5
dimensions = [3,3,3,6]
tpx = rhpx.TaperedHyperX(len(dimensions), dimensions, 1, 1, taper)
# next generate the traffic matrices

