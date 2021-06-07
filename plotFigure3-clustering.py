# -*- coding: utf-8 -*-

import glob, os, sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.style as style
import seaborn as sns
import pandas as pd
from collections import OrderedDict
from matplotlib.colors import LogNorm

fileDirectory = "/home/karpouzi/Research/Eukaryote-mountdir/"

sns.set_context("paper")

wList = [0.0,4.0]
colorList=["b","m","y","r"]
filePaths = []

for w in wList:
    filePath = fileDirectory + "REPS_AGRE_T_4000.0_dtp_0.1_n_40.0_a0_0.2_d0_0.0_ksi_1.2_sar_0.25_a_0.0_w_" + str(w) + "_Tag_0.1_Tab_50.0_Tr_5.0_Td_50.0_d_1.0_dtsave_1.0.csv"
    filePaths.append(filePath)

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))

#space foor w=0 before transition
vec0=np.arange(0,220,30)
vec2=np.arange(950,1150,30)
# plot for w=0
dataFrame = pd.read_csv(filePaths[0], sep=" ", header=0, usecols = ["t","maxSize","nFrag"], nrows=8000)
dataFrame = dataFrame[dataFrame['t']<1100]
dataFrame['t'] = np.around(dataFrame['t'], decimals=0)

df0=dataFrame.loc[dataFrame["t"].isin(vec0)]
df00=df0.groupby(['t']).mean(numeric_only=True)
df1=dataFrame.loc[(dataFrame["t"]>220)&(dataFrame["t"]<950)]
df11=df1.groupby(['t']).mean(numeric_only=True)
df2=dataFrame.loc[dataFrame["t"].isin(vec2)]
df22=df2.groupby(['t']).mean(numeric_only=True)

df=dataFrame.groupby(['t']).mean(numeric_only=True)

sns.set_palette("colorblind")
sns.set_color_codes("colorblind")
sns.scatterplot(x="t",y="nFrag",size="maxSize", color='w', sizes=(df['maxSize'].min()*200, df['maxSize'].max()*200), ax=axs[0,1], data=df)
sns.scatterplot(x="t",y="nFrag",size="maxSize", color='b', sizes=(df00['maxSize'].min()*200, df00['maxSize'].max()*200), ax=axs[0,1], data=df00,legend=False)
sns.scatterplot(x="t",y="nFrag",size="maxSize", color='b', sizes=(df11['maxSize'].min()*200, df11['maxSize'].max()*200), ax=axs[0,1], data=df11,legend=False, edgecolor='face')
sns.scatterplot(x="t",y="nFrag",size="maxSize", color='b', sizes=(df22['maxSize'].min()*200, df22['maxSize'].max()*200), ax=axs[0,1], data=df22,legend=False)


#space foor w=4.0
vec0=np.arange(0,1150,30)
# plot for w=0
dataFrame = pd.read_csv(filePaths[1], sep=" ", header=0, usecols = ["t","maxSize","nFrag"], nrows=8000)
dataFrame = dataFrame[dataFrame['t']<1100]
dataFrame['t'] = np.around(dataFrame['t'], decimals=0)

df0=dataFrame.loc[dataFrame["t"].isin(vec0)]
df00=df0.groupby(['t']).mean(numeric_only=True)
sns.scatterplot(x="t",y="nFrag",size="maxSize", color='r', sizes=(df00['maxSize'].min()*200, df00['maxSize'].max()*200), ax=axs[1,1], data=df00,legend=False)

axs[0,1].annotate(r"$\omega=0.0$", xy=(1010,140))
axs[1,1].annotate(r"$\omega=4.0$", xy=(1010,140))

axs[0,1].set_title(r"$\bf{b}$   Fragmentation measures over time")

axs[0,1].set_ylim(-10,160)
axs[1,1].set_ylim(-10,160)
axs[0,1].legend(title="Maximum size")
axs[0,1].set_xlabel("")
axs[1,1].set_xlabel("Time")
axs[0,1].set_ylabel("Number of fragments")
axs[1,1].set_ylabel("Number of fragments")

dataFrameTheo = pd.read_csv("/home/karpouzi/Research/Eukaryote-mountdir/experimentA0W-220421.csv")
wList = [0.0,1.0,2.0,4.0]

for ix, w in enumerate(wList):
    sns.lineplot(x="N",y="maxSize", ax=axs[0,0], color = colorList[ix], linewidth=2.0 ,label = str(w), data=dataFrameTheo[dataFrameTheo['w']==w] )
    sns.lineplot(x="N",y="nFrag", ax=axs[1,0], color = colorList[ix], linewidth=2.0 ,data=dataFrameTheo[dataFrameTheo['w']==w] )

axs[0,0].set_title(r"$\bf{a}$   Fragmentation measures as a function of the fraction of natural land")
axs[0,0].set_xlabel("")
axs[1,0].set_xlabel("Fraction of natural land")
axs[0,0].set_ylabel("Size of the largest fragment")
axs[1,0].set_ylabel("Number of fragments")
axs[0,0].invert_xaxis()
axs[1,0].invert_xaxis()
axs[0,0].legend(title = r'$\omega$')
plt.tight_layout()


plt.savefig('Figure3-clustering.pdf', format='pdf', dpi = 600, bbox_inches='tight')

plt.show()



# filePathextra = "/home/karpouzi/Research/Eukaryote-mountdir/REPS_AGRE_T_100.0_dtp_0.1_n_40.0_a0_0.35_d0_0.0_ksi_1.2_sar_0.25_a_0.0_w_0.0_Tag_0.1_Tab_50.0_Tr_5.0_Td_50.0_d_1.0_dtsave_0.01.csv"
# dfextra = pd.read_csv(filePathextra, sep=" ", header=0, usecols = ["t","N","maxSize","nFrag"])
# dfextra['N'] = np.around(dfextra['N'], decimals=2)


#
# sns.lineplot(x="N",y="maxSize", ax=axs[0], color = "tab:purple", data=dfextra )
# sns.lineplot(x="N",y="nFrag", ax=axs[1], color = "tab:purple", data=dfextra )

# for ix, file in enumerate(filePaths):
#     print("begin import file " + str(ix) )
#     dataFrame = pd.read_csv(file, sep=" ", header=0, usecols = ["t","N","maxSize","nFrag"], nrows = 12000)
#     # dataFrame=dataFrame[(dataFrame['t']>1500) & (dataFrame['t']<3000)]
#     dataFrame['t'] = np.around(dataFrame['t'], decimals=0)
#     sns.lineplot(x="t",y="N", ax=axs[0], label = wList[ix], data=dataFrame)
#     sns.lineplot(x="t",y="maxSize", ax=axs[1], label = wList[ix], data=dataFrame)
#     # dataFrame['N'] = np.around(dataFrame['N'], decimals=2)
#     print("import done")
#
#     # dataFrame=dataFrame.iloc[:6000,:]
#     print(dataFrame)

    # # plot of connectance over time
    # sns.lineplot(x="t",y="maxSize", ax=axs[0,0], label = wList[ix], data=dataFrame)
    # sns.lineplot(x="t",y="nFrag", ax=axs[1,0], label = wList[ix], data=dataFrame)
    # plot of connectance versus natural land
    # sns.lineplot(x="N",y="maxSize", ax=axs[0], color = colorList[ix], data=dataFrameTheo[dataFrameTheo['w']==wList[ix]] )
    # sns.lineplot(x="N",y="maxSize", ax=axs[0], color = colorList[ix], label = wList[ix], data=dataFrame)
    # sns.lineplot(x="N",y="nFrag", ax=axs[1], color = colorList[ix], data=dataFrameTheo[dataFrameTheo['w']==wList[ix]] )
    # sns.lineplot(x="N",y="nFrag", ax=axs[1], color = colorList[ix], label = wList[ix], data=dataFrame)


# aList = [0.25,0.05,0.1,0.15]
# colorList=["b","m","y","r"]
# filePaths = []
#
# for a in aList:
#     filePath = fileDirectory + "REPS_AGRE_T_3000.0_dtp_0.1_n_40.0_a0_0.25_d0_0.0_ksi_1.2_sar_0.25_a_"+str(a)+"_w_0.0_Tag_0.1_Tab_50.0_Tr_5.0_Td_50.0_d_1.0_dtsave_1.0.csv"
#     filePaths.append(filePath)
#
# fig, axs = plt.subplots(nrows=1, ncols=1)
#
# for ix, file in enumerate(filePaths):
#     dataFrame = pd.read_csv(file, sep=" ", header=0, usecols = ["t","N"])
#     dataFrame['t'] = np.around(dataFrame['t'], decimals=0)
#     sns.lineplot(x="t",y="N", ax=axs, label = aList[ix], data=dataFrame)
#
