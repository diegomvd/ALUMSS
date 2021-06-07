# -*- coding: utf-8 -*-

import glob, os, sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.style as style
import seaborn as sns
from collections import OrderedDict
import pandas as pd

pathMu = "/home/karpouzi/Research/Eukaryote-mountdir/gillespieSensitivityResults5/mu.csv"
pathMuStar = "/home/karpouzi/Research/Eukaryote-mountdir/gillespieSensitivityResults5/muStar.csv"
pathSigma = "/home/karpouzi/Research/Eukaryote-mountdir/gillespieSensitivityResults5/sigma.csv"

dfMu = pd.read_csv(pathMu, sep=",")
dfMuStar = pd.read_csv(pathMuStar, sep=",")
dfSigma = pd.read_csv(pathSigma, sep=",")

# normalize by maximum value for mu: this preserves the sign as a valuable inormation on the direction of changes
dfMu.iloc[:,1:] = dfMu.iloc[:,1:].div(dfMu.iloc[:,1:].abs().max(axis=1),axis='index')
# compute z-scores for muStar: this allows to see the relative importance of each factor namely which ones are more important than the mean
dfMuStar.iloc[:,1:] = dfMuStar.iloc[:,1:].sub(dfMuStar.mean(axis=1),axis='index').div(dfMuStar.std(axis=1),axis='index')
# normalize by maximum value for sigma
dfSigma.iloc[:,1:] = dfSigma.iloc[:,1:].div(dfSigma.iloc[:,1:].abs().max(axis=1),axis='index')


# initialize the axis for first figure: heatmaps
fig, axs = plt.subplots(nrows=1, ncols=3)

# diverging palettes for Mu and MuStar
cmap1 = sns.color_palette("icefire", as_cmap=True)
# sequential palette for Sigma
cmap2 = sns.color_palette("rocket", as_cmap=True)

dfMu_i=dfMu.set_index('output')
dfMuStar_i=dfMuStar.set_index('output')
dfSigma_i=dfSigma.set_index('output')

# removing ripley and meanEs, stdEs from figure to simplify visualization
dfMu_i=dfMu_i.loc[["P","N","A0","A1","D","nFrag","maxSize","stdSize"]]
dfMuStar_i=dfMuStar_i.loc[["P","N","A0","A1","D","nFrag","maxSize","stdSize"]]
dfSigma_i=dfSigma_i.loc[["P","N","A0","A1","D","nFrag","maxSize","stdSize"]]

# transpose to have parameters on the y-axis
dfMu_i = dfMu_i.transpose()
dfMuStar_i = dfMuStar_i.transpose()
dfSigma_i = dfSigma_i.transpose()

# rename columns and indexes to have them match the text
dfMu_i2 = dfMu_i.rename( index={'ksi': r'$y_1$','sar': r'$z$','y0': r'$y_0$','a0': r'$a_0$','w': r'$\omega$','Tag': r'$1/\sigma$','Tr': r'$1/\rho_R$','Td': r'$1/\rho_D$','Tab': r'$1/\rho_L$','a': r'$\alpha$'  })
dfMuStar_i2 = dfMuStar_i.rename( index={'ksi': r'$y_1$','sar': r'$z$','y0': r'$y_0$','a0': r'$a_0$','w': r'$\omega$','Tag': r'$1/\sigma$','Tr': r'$1/\rho_R$','Td': r'$1/\rho_D$','Tab': r'$1/\rho_L$','a': r'$\alpha$'  })
dfSigma_i2 = dfSigma_i.rename( index={'ksi': r'$y_1$','sar': r'$z$','y0': r'$y_0$','a0': r'$a_0$','w': r'$\omega$','Tag': r'$1/\sigma$','Tr': r'$1/\rho_R$','Td': r'$1/\rho_D$','Tab': r'$1/\rho_L$','a': r'$\alpha$'  })

sns.heatmap(data=dfMu_i2,cmap=cmap1,linewidth=1.0,linecolor="white",center=0,cbar_kws={'label': r'Normalized $\mu$',"orientation": "horizontal"},ax=axs[0])
sns.heatmap(data=dfMuStar_i2,cmap=cmap1,linewidth=1.0,linecolor="white",center=0,yticklabels=False,cbar_kws={'label': r'z-scores of $\mu^{\star}$',"orientation": "horizontal"},ax=axs[1])
sns.heatmap(data=dfSigma_i2,cmap=cmap2,linewidth=1.0,linecolor="white",yticklabels=False,cbar_kws={'label': r'Normalized $\sigma$',"orientation": "horizontal"},ax=axs[2])

axs[0].set_xlabel('')
axs[1].set_xlabel('')
axs[2].set_xlabel('')

axs[0].set_title(r"$\bf{a}$")
axs[1].set_title(r"$\bf{b}$")
axs[2].set_title(r"$\bf{c}$")

# initialize the axis for first figure: heatmaps without tr, tab
fig, axs = plt.subplots(nrows=1, ncols=3)

dfMu = pd.read_csv(pathMu, sep=",")
dfMuStar = pd.read_csv(pathMuStar, sep=",")
dfSigma = pd.read_csv(pathSigma, sep=",")

dfMu=dfMu.loc[:,["output","ksi","sar","a","a0","y0","Td","Tag","w"]]
dfMuStar=dfMuStar.loc[:,["output","ksi","sar","a","a0","y0","Td","Tag","w"]]
dfSigma=dfSigma.loc[:,["output","ksi","sar","a","a0","y0","Td","Tag","w"]]

# normalize by maximum value for mu: this preserves the sign as a valuable inormation on the direction of changes
dfMu.iloc[:,1:] = dfMu.iloc[:,1:].div(dfMu.iloc[:,1:].abs().max(axis=1),axis='index')
# compute z-scores for muStar: this allows to see the relative importance of each factor namely which ones are more important than the mean
dfMuStar.iloc[:,1:] = dfMuStar.iloc[:,1:].sub(dfMuStar.mean(axis=1),axis='index').div(dfMuStar.std(axis=1),axis='index')
# normalize by maximum value for sigma
dfSigma.iloc[:,1:] = dfSigma.iloc[:,1:].div(dfSigma.iloc[:,1:].abs().max(axis=1),axis='index')

dfMu_i=dfMu.set_index('output')
dfMuStar_i=dfMuStar.set_index('output')
dfSigma_i=dfSigma.set_index('output')

# removing ripley and meanEs, stdEs from figure to simplify visualization
dfMu_i=dfMu_i.loc[["P","N","A0","A1","D","nFrag","maxSize","stdSize"]]
dfMuStar_i=dfMuStar_i.loc[["P","N","A0","A1","D","nFrag","maxSize","stdSize"]]
dfSigma_i=dfSigma_i.loc[["P","N","A0","A1","D","nFrag","maxSize","stdSize"]]

# transpose to have parameters on the y-axis
dfMu_i = dfMu_i.transpose()
dfMuStar_i = dfMuStar_i.transpose()
dfSigma_i = dfSigma_i.transpose()

dfMu_i2 = dfMu_i.rename( index={'ksi': r'$y_1$','sar': r'$z$','y0': r'$y_0$','a0': r'$a_0$','w': r'$\omega$','Tag': r'$1/\sigma$','Td': r'$1/\rho_D$','a': r'$\alpha$'  })
dfMuStar_i2 = dfMuStar_i.rename( index={'ksi': r'$y_1$','sar': r'$z$','y0': r'$y_0$','a0': r'$a_0$','w': r'$\omega$','Tag': r'$1/\sigma$','Td': r'$1/\rho_D$','a': r'$\alpha$'  })
dfSigma_i2 = dfSigma_i.rename( index={'ksi': r'$y_1$','sar': r'$z$','y0': r'$y_0$','a0': r'$a_0$','w': r'$\omega$','Tag': r'$1/\sigma$','Td': r'$1/\rho_D$','a': r'$\alpha$'  })

sns.heatmap(data=dfMu_i,cmap=cmap1,linewidth=1.0,linecolor="white",center=0,cbar_kws={'label': r'Normalized $\mu$',"orientation": "horizontal"},ax=axs[0])
sns.heatmap(data=dfMuStar_i,cmap=cmap1,linewidth=1.0,linecolor="white",center=0,yticklabels=False,cbar_kws={'label': r'z-scores of $\mu^{\star}$',"orientation": "horizontal"},ax=axs[1])
sns.heatmap(data=dfSigma_i,cmap=cmap2,linewidth=1.0,linecolor="white",yticklabels=False,cbar_kws={'label': r'Normalized $\sigma$',"orientation": "horizontal"},ax=axs[2])


# initialize the axis for second figure: scatter categorical plot
# fig, axs = plt.subplots(nrows=1, ncols=2)
#
# sns.set_palette("colorblind")
# sns.set_color_codes("colorblind")

# sns.catplot(x="day", y="total_bill", kind="swarm", data=tips, ax=axs[0])

plt.show()
