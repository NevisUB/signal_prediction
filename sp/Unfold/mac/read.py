import ROOT

from ROOT import sp

sp.LoadSP()


tf = ROOT.TFile("bin/unfold_data.root","READ")

e = tf.Get("Parameters/Energy")
e.dump()
print e._a

e = tf.Get("Responses/response_0_0_NuMomT_Energy")
e.dump()
