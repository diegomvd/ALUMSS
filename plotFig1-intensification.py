# -*- coding: utf-8 -*-

import glob, os, sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.style as style
import seaborn as sns
import pandas as pd
from collections import OrderedDict
from scipy import interpolate
from matplotlib.gridspec import GridSpec

# path for the data files
path1 = "/home/karpouzi/Research/Eukaryote-mountdir/REPS_AGRE_T_4000.0_dtp_0.1_n_40.0_a0_0.1_d0_0.0_ksi_1.2_sar_0.25_a_0.0_w_0.0_Tag_0.1_Tab_50.0_Tr_5.0_Td_50.0_d_1.0_dtsave_1.0.csv"
path2 = "/home/karpouzi/Research/Eukaryote-mountdir/REPS_AGRE_T_4000.0_dtp_0.1_n_40.0_a0_0.1_d0_0.0_ksi_1.2_sar_0.25_a_0.2_w_0.0_Tag_0.1_Tab_50.0_Tr_5.0_Td_50.0_d_1.0_dtsave_1.0.csv"

replicationData1 = pd.read_csv(path1, sep=" ", header=0)
replicationData2 = pd.read_csv(path2, sep=" ", header=0)

# rounding time for cleaner plots
replicationData1['t'] = np.around(replicationData1['t'], decimals=0)
replicationData2['t'] = np.around(replicationData2['t'], decimals=0)

print(replicationData1)
# replicationData=replicationData.iloc[:3000,:]

sns.set_palette("colorblind")
sns.set_color_codes("colorblind")

# initialize the axis
fig = plt.figure(constrained_layout=True, figsize=(9.5,7.5))
gs = GridSpec(4, 2, figure=fig)
ax00 = fig.add_subplot(gs[0, 0])
ax10 = fig.add_subplot(gs[1, 0])
ax20 = fig.add_subplot(gs[2, 0])
ax30 = fig.add_subplot(gs[3, 0])
ax01 = fig.add_subplot(gs[0, 1])
ax11 = fig.add_subplot(gs[1, 1])
ax231 = fig.add_subplot(gs[2:, 1])


# plot population dynamics
sns.lineplot(x="t",y="P", color='b', ax=ax00, linewidth=3.0, data=replicationData1[replicationData1["t"]<2600])
land_types1=['N','D','A0']
land_colors1=["g", "r", "y"]
for ix, land_type in enumerate(land_types1):
    sns.lineplot(x='t',y=land_type,color=land_colors1[ix],label=land_types1,ax=ax10, linewidth=3.0, data=replicationData1[replicationData1["t"]<2600])
ax00.set_xlabel("")
ax10.set_xlabel("Time")
ax00.set_ylabel("Population density")
ax10.set_ylabel("Fraction of land")
ax00.set_title(r"$\bf{a}$" + "      No intensification "+ r"$(\alpha=0)$")
# ax10.axhline(0.6,linewidth=4, color='k', alpha=0.4)
# ax10.legend(['N','D','A0'])
ax10.legend([r'$N$',r'$D$',r'$A_L$'])

# fig, axs = plt.subplots(nrows=2, ncols=1)
sns.lineplot(x="t",y="P", color='b', ax=ax20, linewidth=3.0, data=replicationData2)
# plot land dynamics
land_types2=['N','D','A0','A1']
land_colors2=["g", "r", "y", "m"]
for ix, land_type in enumerate(land_types2):
    sns.lineplot(x='t',y=land_type,color=land_colors2[ix],label=land_type,ax=ax30, linewidth=3.0, data=replicationData2)
ax20.set_xlabel("")
ax30.set_xlabel("Time")
ax20.set_ylabel("Population density")
ax30.set_ylabel("Fraction of land")
ax20.set_title(r"$\bf{c}$" + "      Low intensification "+ r"$(\alpha=0.2)$")
# ax10.axhline(0.6,linewidth=4, color='k', alpha=0.4)
# ax20.legend(['N','D','A0','A1'])
ax30.legend([r'$N$',r'$D$',r'$A_L$',r'$A_H$'])

###############################################################################

path = "/home/karpouzi/Research/Eukaryote-mountdir/experimentTagA-260421.csv"

def dfStringToList(string_):
    string_ = string_.replace(',','","')
    string_ = string_.replace('[','["')
    string_ = string_.replace(']','"]')
    return string_

# read csv file
samplingData = pd.read_csv(path, sep=",", header=0)
# convert strings to lists in the needed columns
samplingData["N"] = samplingData["N"].apply(dfStringToList)
samplingData["N"] = samplingData["N"].apply(eval)
samplingData["P"] = samplingData["P"].apply(dfStringToList)
samplingData["P"] = samplingData["P"].apply(eval)
samplingData["nMin"] = samplingData["nMin"].apply(dfStringToList)
samplingData["nMin"] = samplingData["nMin"].apply(eval)
samplingData["pMin"] = samplingData["pMin"].apply(dfStringToList)
samplingData["pMin"] = samplingData["pMin"].apply(eval)
samplingData["nMax"] = samplingData["nMax"].apply(dfStringToList)
samplingData["nMax"] = samplingData["nMax"].apply(eval)
samplingData["pMax"] = samplingData["pMax"].apply(dfStringToList)
samplingData["pMax"] = samplingData["pMax"].apply(eval)

# plot for a=0
# fig, axs = plt.subplots(nrows=2, ncols=1)

# plot the SS for P
plotData = samplingData
plotData = plotData.explode("P",ignore_index=True)
plotData["P"] = plotData["P"].astype('float')
plotData=plotData[(plotData['a']==0.0) & (plotData['Tag']>=0.25)]
plotData["Tag"]=1/plotData["Tag"]

sns.lineplot(x="Tag",y="P",ax=ax01,color='b',label='Steady State',linewidth=3,data=plotData)

# # plot the LC for P
plotData = samplingData  #[(samplingData['a']==0.0) & (samplingData['Tag']<0.25) & (samplingData['Tag']>0.025)]
plotData.loc[(plotData["Tag"]<=0.25),'P'] = plotData.loc[(plotData["Tag"]<=0.25),'pMin']
plotData = plotData.explode("P",ignore_index=True)
plotData["P"] = plotData["P"].astype('float')
plotData=plotData[(plotData['a']==0.0) & (plotData['Tag']>0.025) & (plotData['Tag']<=0.25) ]
plotData["Tag"]=1/plotData["Tag"]

sns.lineplot(x="Tag",y="P",ax=ax01,color='m',label='Cycles',linewidth=3,data=plotData)

plotData = samplingData
plotData.loc[(plotData["Tag"]<=0.25),'P'] = plotData.loc[(plotData["Tag"]<=0.25),'pMax']
plotData = plotData.explode("P",ignore_index=True)
plotData["P"] = plotData["P"].astype('float')
plotData=plotData[(plotData['a']==0.0) & (plotData['Tag']>0.025) & (plotData['Tag']<=0.25) ]
plotData["Tag"]=1/plotData["Tag"]
sns.lineplot(x="Tag",y="P",ax=ax01,color='m',linewidth=3,data=plotData)

# plot the SS for N
plotData = samplingData[(samplingData['a']==0.0) & (samplingData['Tag']>=0.25)]
plotData = plotData.explode("N",ignore_index=True)
plotData["N"] = plotData["N"].astype('float')
plotData["Tag"]=1/plotData["Tag"]
sns.lineplot(x="Tag",y="N",ax=ax11,color='b',linewidth=3,data=plotData)
# # plot the LC for N
plotData = samplingData[(samplingData['a']==0.0) & (samplingData['Tag']<0.25) & (samplingData['Tag']>0.025)]
# plotData.loc[(plotData["Tag"]<0.25),'nMin'] = plotData.loc[(plotData["Tag"]<0.25),'N']
plotData = plotData.explode("nMin",ignore_index=True)
plotData["nMin"] = plotData["nMin"].astype('float')
plotData["Tag"]=1/plotData["Tag"]
sns.lineplot(x="Tag",y="nMin",ax=ax11,color='m',linewidth=3,data=plotData)
plotData = samplingData[(samplingData['a']==0.0) & (samplingData['Tag']<0.25) & (samplingData['Tag']>0.025)]
# plotData.loc[(plotData["Tag"]<=0.25),'N'] = plotData.loc[(plotData["Tag"]<=0.25),'nMax']
plotData = plotData.explode("nMax",ignore_index=True)
plotData["nMax"] = plotData["nMax"].astype('float')
plotData["Tag"]=1/plotData["Tag"]
sns.lineplot(x="Tag",y="nMax",ax=ax11,color='m',linewidth=3,data=plotData)

ax01.legend(title= 'Equilibria')
# ax01.annotate(r'$\alpha=0$',xy=(1,320))

ax01.set_xscale("log")
ax11.set_xscale("log")

ax01.set_xlabel('')
ax11.set_xlabel('Responsiveness to resource demand '+r"$\sigma$")
ax11.set_ylabel('Natural land fraction')
ax01.set_ylabel('Population density')

ax01.set_title(r"$\bf{b}$"+"        Bifurcation diagrams for "+r"$\alpha=0$")
# ax11.axhline(0.59,linewidth=4, color='k', alpha=0.4)

###############################################################################

###############################################################################
# 2d bifurcation diagram

path = "/home/karpouzi/Research/Eukaryote-mountdir/experimentTagA-290421-2Dbif.csv"

# read csv file
df = pd.read_csv(path, sep=",", header=0)
df = df.loc[(df["a"]<1.0) & (df["a"]>0.1)]
alist=df['a'].unique()

# identify the threshold in Tag for each a
# get the points at collapse
df1 = df[(df['P']<10.0) & (df['N']==0.0)]

tag_th1=np.zeros(len(alist))
for ix,a in enumerate(alist):
    df3 = df1.loc[df1['a']==a]
    tag_th1[ix] = df3['Tag'].max()

# get the points at sustainability
df2 = df[(df['P']<10.0) & (df['N']==0.0)]
tag_th2=np.zeros(len(alist))
for ix,a in enumerate(alist):
    df3 = df2.loc[df2['a']==a]
    tag_th2[ix] = df3['Tag'].min()

tag_th = tag_th1*0.5 + tag_th2*0.5
# adding the point for a=0.1 from other simulations
tag_th = np.append(tag_th,0.15)
alist = np.append(alist,0.1)
# adding the hopf critical point for a=0.0
tag_th = np.append(tag_th,0.2)
alist = np.append(alist,0.0)
# adding extra point for correct fill in sustainable region
tag_th = np.append(tag_th,1.0)
alist = np.append(alist,0.0)

tag_th=1/tag_th

icefire=sns.color_palette("icefire",2)
sns.set_palette(icefire)
# sns.set(palette='icefire',color_codes=True)
# plot the scatter points
# fig, axs = plt.subplots(nrows=1, ncols=1)
sns.scatterplot(x=tag_th,y=alist,color='k',ax=ax231)
sns.lineplot(x=[1/0.0178,1/0.19],y=[0.0,0.0],color='m',linewidth=4,ax=ax231)


f = interpolate.interp1d(tag_th, alist)
tag_new = np.arange(np.min(tag_th), max(tag_th), 0.01)
a_new = f(tag_new)

ax231.fill_between(tag_new,a_new,color='r',alpha=0.7)
# ax231.fill_betweenx(a_new,tag_new,0.01,color='r',alpha=0.3)
# ax231.fill_betweenx([0.7,0.8,0.9],[0.018,0.018,0.018],0.01,color='r',alpha=0.3)
ax231.fill_between(tag_new,a_new,0.9,color='b',alpha=0.7)

ax231.annotate("Sustainable steady state",xy=(2.0,0.5))
ax231.annotate("Irreversible \n  collapse",xy=(16,0.15))
ax231.annotate("Reversible collapses \n         (cycles)",xy=(9.0,0.0),xytext=(1.1,0.1),arrowprops=dict(facecolor='black', shrink=1.0, width=2.0, headwidth=6.5, alpha=0.7))
ax231.set_xscale('log')
ax231.set_xlabel('Responsiveness to resource demand '+r"$\sigma$")
ax231.set_ylabel('Preference for intensification '+r"$\alpha$")
ax231.set_title(r"$\bf{d}$" + "      2D Bifurcation diagram")
#
# plot for a>0
# fig, axs = plt.subplots(nrows=2, ncols=1)
#
# palette=sns.color_palette("flare")
#
# plotData = samplingData[samplingData['a']>0.0]
# plotData = plotData.explode("P",ignore_index=True)
# plotData["P"] = plotData["P"].astype('float')
# plotData["Tag"]=1/plotData["Tag"]
#
# sns.lineplot(x="Tag",y="P",ax=axs[0],linewidth=2,color=palette[2],label='1/4',data=plotData[plotData['a']==0.25])
# sns.lineplot(x="Tag",y="P",ax=axs[0],linewidth=2,color=palette[4],label='1/2',data=plotData[plotData['a']==0.5])
#
# plotData = samplingData[samplingData['a']>0.0]
# plotData = plotData.explode("N",ignore_index=True)
# plotData["N"] = plotData["N"].astype('float')
# plotData["Tag"]=1/plotData["Tag"]
#
# sns.lineplot(x="Tag",y="N",ax=axs[1],color=palette[2],linewidth=2,data=plotData[plotData['a']==0.25])
# sns.lineplot(x="Tag",y="N",ax=axs[1],color=palette[4],linewidth=2,data=plotData[plotData['a']==0.5])
#
# axs[0].legend(title = r"$\alpha$")
#
# axs[0].set_xscale("log")
# axs[1].set_xscale("log")
#
# axs[0].set_xlabel('')
# axs[1].set_xlabel('Responsiveness to resource demand')
# axs[1].set_ylabel('Fraction of natural land')
# axs[0].set_ylabel('Population density')
#
# axs[0].set_title(r"$\bf{d}$" + "      Intensification")
# axs[1].axhline(0.6,linewidth=4, color='k', alpha=0.4)

plt.savefig('Figure1-intensification-revision.pdf', format='pdf', dpi = 1200, bbox_inches='tight')

plt.show()
