from ROOT import sp
sp.LoadSP()

a = sp.SPIO()
a.add_mc_in_file("/Users/vgenty/Desktop/signal_prediction/simplifyTreeOsc/output_osc_mc_detail_1.root")
a.set_mc_tree_name("MiniBooNE_CCQE")

a.initialize()

a.dump_branches()
