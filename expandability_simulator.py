import sys, math
import hyperx as hx
import dragonfly as df

eps_radix = 32

# figure out what the expandability of a fixed radixed network is
# for hyperX

S_2D = int(math.floor((eps_radix ** (1. / 2)) + 1)) # number of switches in each dimension assuming that each dimension has the same number of switches
S_3D = int(math.floor((eps_radix ** (1. / 3)) + 1))
S_4D = int(math.floor((eps_radix ** (1. / 4)) + 1))

D2HyperX = hx.HyperX(2, [7, 6], 1,1)
D3HyperX = hx.HyperX(3, [6, 4, 3], 1,1)
D4HyperX = hx.HyperX(4, [4,3,3,3], 1,1)

num_dfly_group = eps_radix + 1
dfly = df.CanonicalDragonfly(num_dfly_group, 1)

print "EPS radix: {}".format(eps_radix)
print "2D HyperX: {}".format(D2HyperX.cheeger_constant())
print "3D HyperX: {}".format(D3HyperX.cheeger_constant())
print "4D HyperX: {}".format(D4HyperX.cheeger_constant())
#print "Dragonfly: {}".format(dfly.cheeger_constant())


eps_radix = 64

# figure out what the expandability of a fixed radixed network is
# for hyperX

S_2D = int(math.floor((eps_radix ** (1. / 2)) + 1)) # number of switches in each dimension assuming that each dimension has the same number of switches
S_3D = int(math.floor((eps_radix ** (1. / 3)) + 1))
S_4D = int(math.floor((eps_radix ** (1. / 4)) + 1))

D2HyperX = hx.HyperX(2, [9, 9], 1,1)
D3HyperX = hx.HyperX(3, [6, 4, 3], 1,1)
D4HyperX = hx.HyperX(4, [4,3,3,3], 1,1)

num_dfly_group = eps_radix + 1
dfly = df.CanonicalDragonfly(num_dfly_group, 1)

print "EPS radix: {}".format(eps_radix)
print "2D HyperX: {}".format(D2HyperX.cheeger_constant())
print "3D HyperX: {}".format(D3HyperX.cheeger_constant())
print "4D HyperX: {}".format(D4HyperX.cheeger_constant())
print "Dragonfly: {}".format(dfly.cheeger_constant())