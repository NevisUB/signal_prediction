import sys, os

import ROOT
from ROOT import sp
sp.LoadCombined()

import matplotlib
from matplotlib import colors as mcolors
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica Neue']

import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

import root_numpy as rn

df = pd.DataFrame(rn.root2array([sys.argv[1]],treename='MiniBooNE_CCQE'))


df['Eqe'] = df.apply(lambda x : sp.CombinedFit_EnuQE_ryan(x['Energy'],x['CosTheta'],3),axis=1)


print "...Computing Pi0..."
df['Pi0'] = df.apply(lambda x : sp.Pi0Details(x['NFSP'],
                                              x['FSPType'],
                                              x['VertexX'],x['VertexY'],x['VertexZ'],
                                              x['MomX'],x['MomY'],x['MomZ'],x['MomT']),axis=1)

print "...Computing Backgrounds..."
df['Stacked'] = df.apply(lambda x : sp.StackHistoBkgd(0,bool(x['Pi0']),x['NUANCEChan'],x['NuType'],x['NuParentID']),axis=1)

print df['Stacked']

bins=np.array([200,300,375,475,550,675,800,950,1100,1300,1500,3000])

data_v   = []
weight_v = []
hist_v   = []

#bkgd_v = ['kBKGD_INVALID','kBKGD_DIRT','kBKGD_PI0','kBKGD_DELTA','kBKGD_NUEPIP','kBKGD_NUEKP','kBKGD_NUEK0','kBKGD_OTHER']

bkgd_v  = ['kBKGD_PI0',
           'kBKGD_DELTA',
           'kBKGD_NUEPIP',
           'kBKGD_NUEKP',
           'kBKGD_NUEK0',
           'kBKGD_OTHER']

names_v = [r"$\pi^0$ misid",
           r"$\Delta\rightarrow$N$\gamma$",
           r"$\nu_e$ from $\mu^+$",
           r"$\nu_e$ from K$^+$",
           r"$\nu_e$ from K$^0$",
           "Other"]

color_v = ['#CF5E60',
           '#DDB989',
           '#B4CDC7',
           '#85C1A5',
           '#839E8F',
           '#999999']

for bkgd,bkgd_name in enumerate(bkgd_v):
    bkgd += 2
    
    this_df = df.query("PassOsc==1 & Stacked==@bkgd")
    data    = this_df.Eqe.values*1000.
    #data = this_df.RecoEnuQE.values*1000.
    weight  = this_df.Weight.values*0.157
    hist = np.histogram(data,weights=weight,bins=bins)
    hist = hist[0]
    data_v.append(data)
    weight_v.append(weight)
    hist_v.append(hist)

    print "\tFilled...",bkgd_name
    
   
print "...Drawing..."

fig,ax=plt.subplots(figsize=(10,6))
hist_v=np.array(hist_v)

bin_heights=np.sum(hist_v,axis=0)
dbins=bins[1:]-bins[:-1]

hist_v=np.divide(hist_v,dbins)

bins_cpy = list(bins)
bins_cpy[-1] = 1700.0
a = [bins.copy()[:-1] / 1000.0]*6
b = bins.copy() / 1000.0

hist_v = list(hist_v)
order   = [5,1,0,4,3,2]

#
# dirt
#
dirt_v = np.array([0.114947, 0.0846428, 0.0606367, 0.0492776, 0.030317, 0.0126265, 0.00633687, 0.00760151,0,0,0])
hist_v += [dirt_v]
a = [bins.copy()[:-1] / 1000.0]*7
order   = [5,6,1,0,4,3,2]
color_v += ["#866658"]
names_v += ["Dirt"]


hist_v  = [hist_v[i]  for i in order]
color_v = [color_v[i] for i in order]
names_v = [names_v[i] for i in order]


n,b,p =ax.hist(a,
               bins=b,
               weights=hist_v,
               stacked=True,
               histtype='stepfilled',
               label = names_v,
               color = color_v)

ax.legend(fontsize=20)
ax.set_xlabel(r"$\mathrm{{\bf E}}^{\mathrm{{\bf QE}}}_{\nu}$ (GeV)",
              fontweight='bold',
              fontsize=25,
              horizontalalignment='right', 
              x=1.0)

ax.set_ylabel("Events / MeV",
              fontweight='bold',
              fontsize=25,
              horizontalalignment='right', 
              y=1.0)

ticks=np.array([0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.7])

tick_labels = list(ticks)
tick_labels[-1]=3.0
ax.set_xticks(ticks)
ax.set_xticklabels(tick_labels,fontweight='bold',fontsize=15)

ax.set_xlim(ticks[0],ticks[-1])

ax.tick_params(which = 'minor', direction = 'in', length=5)
ax.tick_params(which = 'major', direction = 'in', length=10)
minor_ticks_x = np.arange(0.2+0.05,1.7,0.05)
ax.set_xticks(minor_ticks_x, minor = True)

ax.set_ylim(0,2.5)
minor_ticks_y = np.arange(0+0.1,2.5,0.1)
ax.set_yticks(minor_ticks_y, minor = True)
ax.set_yticklabels(plt.yticks()[0],fontweight='bold',fontsize=15)

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.tight_layout()
plt.savefig("stacked_background.pdf")

plt.cla()
plt.clf()
plt.close()

