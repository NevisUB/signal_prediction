import ROOT
from ROOT import sp
import numpy as np


sp.LoadSP()

a = sp.SPIO()
a.add_mc_in_file("/Users/vgenty/Desktop/signal_prediction/simplifyTreeOsc/output_osc_mc_detail_1.root")
a.set_mc_tree_name("MiniBooNE_CCQE")

a.initialize()


#
# Evis spectrum
# 
bins_lo_v = np.arange(0,3000,200).astype(np.float64)
bins_hi_v = bins_lo_v+200.0

bins_lo_v = sp.as_vector_double(bins_lo_v)
bins_hi_v = sp.as_vector_double(bins_hi_v)

name = ROOT.std.string("Energy")
a.add_reco_parameter(name,bins_lo_v,bins_hi_v);

#
# E true spectrum
# 
bins_lo_v = np.arange(0,3.0,0.2).astype(np.float64)
bins_hi_v = bins_lo_v+0.2

bins_lo_v = sp.as_vector_double(bins_lo_v)
bins_hi_v = sp.as_vector_double(bins_hi_v)

name = ROOT.std.string("NuMomT")
a.add_true_parameter(name,bins_lo_v,bins_hi_v);

a.init_response_matrix();
a.fill_responses();
a.write_unfold_file();

print "Calling SPIO destructor"
del a
print "bye?"
