import ROOT

from ROOT import sp

sp.LoadSP()


tf = ROOT.TFile("unfold_data.root","READ")

e = tf.Get("Parameters/Energy")
e.dump()


e = tf.Get("Responses/response_0_0_NuMomT_Energy")
e.dump()
