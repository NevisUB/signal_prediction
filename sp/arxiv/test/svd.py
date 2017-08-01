from ROOT import sp
print "...Load unfold..."
sp.LoadSP()
sp.LoadUnfold()
print "...Loaded!..."

import numpy as np

dim=4

a = np.random.rand(dim*dim).reshape(dim,dim).astype(np.float32)
print "Numpy a"
print a
print
tsvd = sp.TikhonovSVD()

print "Mat a"
k = sp.as_mat_float32(a)
print k
print k.rows(),k.cols()
print

tsvd.Initialize(k)

print "as array tsvd.A()"
print sp.as_array_float32(tsvd.A())
print
print
np.set_printoptions(3)
print np.abs(sp.as_array_float32(tsvd.A())-a)
print
print "as array k"
print sp.as_array_float32(k) 
print 
import ROOT

e = ROOT.TH2F("","",4,0,4,2,0,4)

e.Fill(1,1)
e.Fill(2,1)
e.Fill(3,1)

print e.GetNbinsX(),e.GetNbinsY(),e.GetSize()

print e
b = sp.to_mat_eigen(e)
print b
print sp.as_array_float32(b)
e.Draw("COLZ")
raw_input('')
