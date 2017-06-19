import os,sys

import numpy as np
import pandas as pd

import root_numpy as rn

from plot_util import Efficiency1D, StackedBackgrounds

from ROOT import sp
sp.LoadCombined()


FNAME_v=sys.argv[1:]

print "Got input files: ",FNAME_v

print "Reading in events..."
tree_df = pd.DataFrame(rn.root2array(FNAME_v,selection='Weight>0',treename='MiniBooNE_CCQE'))

print "Calculating Eqe..."
tree_df['Eqe']    = tree_df.apply(lambda x : 1000.0*sp.CombinedFit_EnuQE_ryan(x['Energy'],x['CosTheta'],x['NuType']),axis=1)

print "Calculating backgrounds..."
tree_df['BkgdID'] = tree_df.apply(lambda x : sp.CombinedFit_bkgd_type(x['NUANCEChan'],x['NuType']),axis=1)

#
# Signal efficiency
#

# Evis
query1 = "PassOsc==1 & (NuType==3 | NuType==4) & NUANCEChan==1"
query2 = "             (NuType==3 | NuType==4) & NUANCEChan==1"
xlo = 0
xhi = 2000
dx  = 100
Efficiency1D(tree_df,query1,query2,"Energy",xlo,xhi,dx)

query3 = "PassOsc==1"
StackedBackgrounds(tree_df,query3,"Energy",xlo,xhi,dx)

# Eqe
query1 = "PassOsc==1 & (NuType==3 | NuType==4) & NUANCEChan==1"
query2 = "             (NuType==3 | NuType==4) & NUANCEChan==1"
xlo = 0
xhi = 2000
dx  = 100
Efficiency1D(tree_df,query1,query2,"Eqe",xlo,xhi,dx)
