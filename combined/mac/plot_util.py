import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=20
matplotlib.rcParams['font.family']='Times New Roman'

import numpy as np
import pandas as pd

def Efficiency1D(df,query1,query2,var,xlo,xhi,dx,useweight=True):

    df_1 = df.query(query1)
    df_2 = df.query(query2)

    fig,ax=plt.subplots(figsize=(10,6))

    data1   = df_1[var].values
    data2   = df_2[var].values

    weight1 = weight2 = np.ones(data1.size)

    if useweight:
        weight1 = df_1['Weight'].values
        weight2 = df_2['Weight'].values

    bins=np.arange(xlo,xhi+dx,dx)
    thist = ax.hist([data1,data2],
                    weights=[weight1,weight2],
                    bins=bins,
                    stacked=True,
                    histtype='stepfilled',
                    label=['Selected','All'])

    ax.legend(loc='best')
    ax.set_xlabel("Evis",fontweight='bold')
    ax.set_xlim(xlo,xhi)
    ax.grid()
    SS="{}_pass_all.pdf".format(var)
    print "Write {}".format(SS)
    plt.savefig(SS,format='pdf')
    plt.cla()
    plt.clf()
    plt.close() 

    fig,ax=plt.subplots(figsize=(10,6))

    bin_w = (thist[1][1]-thist[1][0])
    bin_dq = bin_w / 2.0

    xdata = thist[1][:-1] + bin_dq
    ydata = thist[0][0] / thist[0][1]

    yerr = np.sqrt(ydata*(1-ydata) / thist[0][1])
    yerr = np.nan_to_num(yerr)
    xerr = np.ones(yerr.size) * bin_dq 
    ax.errorbar(xdata,ydata,xerr=xerr,yerr=yerr,fmt='o',lw=2,color='blue')
    
    print "Efficiency: ",str(ydata)
    ax.grid()
    SS="{}_efficiency.pdf".format(var)
    print "Write {}".format(SS)
    plt.savefig(SS,format='pdf')
    plt.cla()
    plt.clf()
    plt.close() 

    
def StackedBackgrounds(df,query,var,xlo,xhi,dx):
    bins=np.arange(xlo,xhi+dx,dx)
    
    data_v = []
    weight_v = []

    bkgd_v = ['BINVALID','BCCQE','BCCPIP','BNCPI0','BCHPI0','BDELTA']
    color_v = ['red','green','magenta','yellow','blue','cyan']

    for bkgd_id,bkgd_name in enumerate(bkgd_v):

        bkg_query = query + " & BkgdID==@bkgd_id"
        this_df  = df.query(bkg_query)

        data_v.append(this_df[var].values)
        weight_v.append(this_df[var].values)

    fig,ax=plt.subplots(figsize=(10,6))

    thist = ax.hist(data_v,
                    weights = weight_v,
                    color   = color_v,
                    bins    = bins,
                    stacked = True,
                    histtype= 'stepfilled',
                    label   = bkgd_v)

    ax.legend()
    ax.set_xlabel("Evis",fontweight='bold')
    ax.set_xlim(0,2000)
    ax.grid()
    
    SS="{}_background.pdf".format(var)
    print "Write {}".format(SS)
    plt.savefig(SS,format='pdf')
    plt.cla()
    plt.clf()
    plt.close() 
    
