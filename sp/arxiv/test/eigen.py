from ROOT import sp
sp.LoadAlgorithm()

import numpy as np

a = np.random.rand(6).reshape(3,2).astype(np.float32)

print a

k = sp.as_mat_float32(a);

print k
print k.rows()
print k.cols()

print "["
for i in xrange(k.rows()):
    SS = "["
    for j in xrange(k.cols()):
        SS += str(k(i,j)) + ","
    print SS + "]"
print "]"

print
l = sp.as_array_float32(k)

print l

