# -*- coding: utf-8 -*-

import glob, os, sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.style as style
import seaborn as sns
import pandas as pd
from collections import OrderedDict

# path for the data files
path1 = "/home/karpouzi/Research/Eukaryote-mountdir/REPS_AGRE_T_4000.0_dtp_0.1_n_40.0_a0_0.1_d0_0.0_ksi_1.2_sar_0.25_a_0.0_w_0.0_Tag_0.1_Tab_50.0_Tr_5.0_Td_50.0_d_1.0_dtsave_1.0.csv"
path2 = "/home/karpouzi/Research/Eukaryote-mountdir/REPS_AGRE_T_100.0_dtp_0.1_n_40.0_a0_0.35_d0_0.0_ksi_1.2_sar_0.25_a_0.0_w_0.0_Tag_0.1_Tab_50.0_Tr_5.0_Td_50.0_d_1.0_dtsave_0.01.csv"

sns.set_context("paper")

replicationData1 = pd.read_csv(path1, sep=" ", header=0)
replicationData1 = replicationData1.loc[replicationData1['t']<=825]
replicationData1['conversionProb'] = replicationData1['P']-replicationData1['Y']
replicationData1['binN'] = np.around(replicationData1['N'], decimals=2)
replicationData1['binSize'] = np.around(replicationData1['maxSize'], decimals=2)
replicationData1['t'] = np.around(replicationData1['t'], decimals=0)

replicationData2 = pd.read_csv(path2, sep=" ", header=0)
# replicationData2 = replicationData2.loc[replicationData2['t']<=100]
replicationData2['conversionProb'] = replicationData2['P']-replicationData2['Y']
replicationData2['binN'] = np.around(replicationData2['N'], decimals=2)
replicationData2['binSize'] = np.around(replicationData2['maxSize'], decimals=2)
replicationData2['t'] = np.around(replicationData2['t'], decimals=0)

# filters for accurate representation of resources and expansion propensity
replicationData2PostPerco = replicationData2.loc[(replicationData2['N']<0.6)]
# this filters the sharp reduction in agro propensity when N=0 and the following state at N<0.1 with propensity=0
replicationData2Agro1 = replicationData2PostPerco.loc[(replicationData2PostPerco['conversionProb']>10) & (replicationData2PostPerco['N']>0.01)]
# now get pre perco to join them
replicationData2PrePerco = replicationData2.loc[(replicationData2['N']>=0.6)]
# now join the two of them
replicationData2Agro1 = pd.concat([replicationData2PrePerco,replicationData2Agro1])
# and now join with data1
replicationDataAgro1 = pd.concat([replicationData1.loc[replicationData1['N']>0.5],replicationData2Agro1])
# now get what was previously filtered out: line at conversionProb = 0
replicationDataAgro2 = replicationData2.loc[(replicationData2['conversionProb']<10) & (replicationData2['N']>0.01) & (replicationData2['N']<0.25)]
replicationDataAgro3 = replicationData2.loc[(replicationData2['N']<=0.01)]
replicationDataAgro3.sort_values(by=['conversionProb'])
# print(replicationDataAgro3)

frames = [replicationData1,replicationData2]
replicationData = pd.concat(frames)

# print(replicationData)

# replicationData=replicationData.iloc[:3000,:]

# # initialize the axis
# fig, axs = plt.subplots(nrows=2, ncols=1)
#
# # plot population dynamics
# sns.lineplot(x="t",y="P", ax=axs[0], data=replicationData)
# print("yaesta")
#
# # plot land dynamics
# land_types=["N","D","A0","A1"]
# land_colors=["tab:green", "tab:red", "tab:orange", "tab:purple"]
# for ix, land_type in enumerate(land_types):
#     sns.lineplot(x='t',y=land_type,color=land_colors[ix],label=land_type,ax=axs[1], data=replicationData)
#     # sns.lineplot(x='t',y=land_type,color=land_colors[ix],label=land_type,ax=axs[1], data=replicationData2)
#
# # setting labels and plotting legend
# axs[1].set_xlabel("Time")
# axs[0].set_ylabel("Population density")
# axs[1].set_ylabel("Fraction of land")
# plt.legend()

###############################################################################
# 2- plotting the fragmentation and ecosystem service metrics

# initialize the axis
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(7.5,4.5))

sns.set_palette("deep")
sns.set_color_codes("deep")

# plot the maximum natural fragment size against the fraction of natural land
# sns.scatterplot(x="N",y="maxSize",color="tab:blue", ax=axs[0,0],data=replicationData)
sns.lineplot(x='binN',y="maxSize",color='b',ax=axs[0,0],linewidth=3.0, data=replicationData)

# plot the number of natural fragments against the fraction of natural land
# sns.scatterplot(x="N", y="nFrag",color="tab:blue",ax=axs[1,0], data = replicationData)
sns.lineplot(x='binN', y="nFrag",color='b',ax=axs[1,0],linewidth=3.0, data = replicationData)

# plot the mean ES provision against the max fragment size
# sns.scatterplot(x="maxSize", y="meanES", color="tab:green", ax=axs[0,1], data=replicationData)
sns.lineplot(x='binSize', y="meanES", color="g",linewidth=3.0, ax=axs[0,1], data=replicationData)

# plot the var ES provision against the max fragment size
# sns.scatterplot(x="maxSize", y="stdES", color="tab:green", ax=axs[1,1], data=replicationData)
sns.lineplot(x='binSize', y="stdES", color="g",linewidth=3.0,  ax=axs[1,1], data=replicationData)

# plot the agricultural production against the fraction of natural land
# averageResource = replicationData['Y']#/replicationData['A0']/1600
# sns.scatterplot(x="N",y=averageResource, color="tab:orange",ax=axs[0,2], data=replicationData)
# sns.scatterplot(x="N",y="Y", color="y",ax=axs[0,2], data=replicationData)
sns.lineplot(x='binN',y="Y", color="y",linewidth=3.0,ax=axs[0,2], data=replicationDataAgro1)
sns.lineplot(x='binN',y="Y", color="y",linewidth=3.0,ax=axs[0,2], data=replicationDataAgro2)

# plot the expansion probability per unit time against fraction of natural land
# sns.scatterplot(x="N", y=conversionProb, color="tab:red",ax=axs[1,2], data=replicationData)
# sns.scatterplot(x="N", y=conversionProb, color="r",ax=axs[1,2], data=replicationData)
sns.lineplot(x='binN', y='conversionProb', color="r",linewidth=3.0,ax=axs[1,2], data=replicationDataAgro1)
sns.lineplot(x='binN', y='conversionProb', color="r",linewidth=3.0,ax=axs[1,2], data=replicationDataAgro2)
axs[1,2].plot([0,0],[2,97.2],color='r',linewidth=3.0)
# axs[1,2].plot(replicationDataAgro3['binN'],replicationDataAgro3['conversionProb'],color='r',linewidth=2.0)
# sns.scatterplot(x='binN', y='conversionProb', color='r', ax=axs[1,2], data=replicationDataAgro3)


# setting labels and legend
# axs[1,2].set_xlim(left=1.0, right = 0.0)
# axs[0,0].set_xlim(left=1.0)
# axs[1,0].set_xlim(left=1.0)
# axs[0,2].set_xlim(left=1.0, right = 0.0)

axs[0,0].set_title(r"$\bf{a}$   Fragmentation measures")
axs[0,1].set_title(r"$\bf{b}$   Ecosystem services (ES)")
axs[0,2].set_title(r"$\bf{c}$   Societal feedback")

axs[0,0].invert_xaxis()
axs[1,0].invert_xaxis()
axs[0,2].invert_xaxis()
axs[1,2].invert_xaxis()

axs[0,0].axvline(0.59,linewidth=4, color='k', alpha=0.4)
axs[1,0].axvline(0.59,linewidth=4, color='k', alpha=0.4)
axs[0,2].axvline(0.59,linewidth=4, color='k', alpha=0.4)
axs[1,2].axvline(0.59,linewidth=4, color='k', alpha=0.4)

axs[0,0].annotate('Percolation \n threshold', xy=(0.55, 0.55),xytext=(0.45, 0.74),arrowprops=dict(facecolor='black', shrink=1.1, width=2.0, headwidth=6.5, alpha=0.7))


axs[0,0].set_ylabel("Largest fragment size")
axs[0,0].set_xlabel("")
axs[1,0].set_ylabel("Number of fragments")
axs[1,0].set_xlabel("Fraction of natural land")

axs[0,1].set_ylabel("Mean ES provision")
axs[0,1].set_xlabel("")
axs[1,1].set_ylabel("SD of ES provision")
axs[1,1].set_xlabel("Largest fragment size")

axs[0,2].set_ylabel("Total resource \n production per unit time")
axs[0,2].set_xlabel("")
axs[1,2].set_ylabel("Agricultural \n expansion propensity")
axs[1,2].set_xlabel("Fraction of natural land")

plt.tight_layout()
plt.savefig('Figure2-fragmentation-revision.pdf', format='pdf', dpi = 1200, bbox_inches='tight')

plt.show()
