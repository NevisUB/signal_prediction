from ROOT import sp
sp.LoadCombined()

import numpy as np

b=[i+5 for i in xrange(100)]
a = np.array(b,dtype=np.float32)
print a
vec = sp.as_vector_float32(a)
print
print "vector size",vec.size()
for i in vec:
    print i


