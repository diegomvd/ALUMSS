# -*- coding: utf-8 -*-

import glob, os, sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.style as style
import seaborn as sns
import pandas as pd
from collections import OrderedDict

# path for the data files
path = "/home/karpouzi/Research/Chapter2/gillespie-land-use/experimentPercoMeasure-280521.csv"

df = pd.read_csv(path, sep=",")

# initialize the axis for first figure: radius of gyration and correlation length as a function of N and for several w
fig, axs = plt.subplots(nrows=1, ncols=3)

# # setting the real scale to plot log-log
# df.loc[:,['corrLen']] = df.loc[:,['corrLen']]*100
# df.loc[:,['N']] = df.loc[:,['N']]*10000
print(df.loc[:,['N']])
df.loc[:,['meanES']] = df.loc[:,['meanES']].div(df.loc[:,['corrLen']],axis=0)
print(df.loc[:,['meanES']])


sns.set_palette("colorblind")
sns.scatterplot(x="N",y="corrLen",ax=axs[0],hue="w",data=df)
sns.scatterplot(x="corrLen",y="meanES",ax=axs[1],hue="w",data=df)
sns.scatterplot(x="N",y="meanES",ax=axs[2],hue="w",data=df)
# axs[0].set_yscale('log')
# axs[1].set_yscale('log')
# axs[2].set_yscale('log')
# axs[0].set_xscale('log')
# axs[1].set_xscale('log')
# axs[2].set_xscale('log')

# # mean ecosystem service provision against correlation length and std ES against correlation length for several w
# fig, axs = plt.subplots(nrows=1, ncols=2)
# sns.scatterplot(x="corrLen",y="meanES",ax=axs[0],hue="w",data=df)
# sns.scatterplot(x="corrLen",y="stdES",ax=axs[1],hue="w",data=df)
# # axs[0].set_xscale('log')
# axs[0].set_yscale('log')
# axs[1].set_yscale('log')


plt.show()
