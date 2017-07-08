from ROOT import sp
sp.LoadCombined()

import numpy as np

dim=4

a = np.random.rand(dim*dim).reshape(dim,dim).astype(np.float32)

tsvd = sp.TikhonovSVD()

k = sp.as_mat_float32(a)
print k
print k.rows(),k.cols()

tsvd.Initialize(k)

print a
print
print
print sp.as_array_float32(tsvd.A())
print
print
np.set_printoptions(3)
print np.abs(sp.as_array_float32(tsvd.A())-a)
print
print sp.as_array_float32(k) 

