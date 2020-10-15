# -*- coding: utf-8 -*-

import glob, os, sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.style as style
import seaborn as sns
from collections import OrderedDict

count=0
filespop=[]
filesland=[]
filesevnt=[]
filesclus=[]
filesclux=[]
for file in glob.glob("DATA_14-10-2020/*.dat"):
    if "DATA_POPU" in file:
        filespop.append(file)
        count+=1
    elif "DATA_LAND" in file:
        filesland.append(file)
        count+=1
    elif "DATA_EVNT" in file:
        filesevnt.append(file)
        count+=1
    elif "DATA_CLUSS" in file:
        filesclus.append(file)
        count+=1
    elif "DATA_CLUX" in file:
        filesclux.append(file)
        count+=1
    else:
        print("Error: for instance there is only population and land files")

print("holaaaaa")

# sort filenames so they match
suffixpop=[]
for file in filespop:
    suffixpop.append(file[file.index('_T_'):])

i=0
while i in range(len(filesland)):
    file=filesland[i]
    ix = suffixpop.index(file[file.index('_T_'):])
    if i!=ix:
        filesland[ix], filesland[i] = filesland[i], filesland[ix]
        i=0
    else:
        i+=1

i=0
while i in range(len(filesevnt)):
    file=filesevnt[i]
    ix = suffixpop.index(file[file.index('_T_'):])
    if i!=ix:
        filesevnt[ix], filesevnt[i] = filesevnt[i], filesevnt[ix]
        i=0
    else:
        i+=1

i=0
while i in range(len(filesclus)):
    file=filesclus[i]
    ix = suffixpop.index(file[file.index('_T_'):])
    if i!=ix:
        filesclus[ix], filesclus[i] = filesclus[i], filesclus[ix]
        i=0
    else:
        i+=1

i=0
while i in range(len(filesclux)):
    file=filesclux[i]
    ix = suffixpop.index(file[file.index('_T_'):])
    if i!=ix:
        filesclux[ix], filesclux[i] = filesclux[i], filesclux[ix]
        i=0
    else:
        i+=1
###############################################################################

sns.set_context('paper')
style.use('seaborn-paper')

for ix in range(count):

    ix1=filespop[ix].find("_a_")+3
    ix2=filespop[ix].find("_w_")

    a=filespop[ix][ix1:ix2]

    ix1=filespop[ix].find("_m_")+3
    ix2=filespop[ix].find("_g_")

    m=filespop[ix][ix1:ix2]

    datapop=np.loadtxt(filespop[ix])
    dataland=np.loadtxt(filesland[ix])
    # dataclus=np.loadtxt(filesclus[ix])

    ## figure 1 population, land and production over time
################################################################################
    fig, axs = plt.subplots(2,1,sharex='col')

    axs[0].plot(datapop[:,0],datapop[:,1])

    colors=['tab:green','tab:orange','tab:red']
    land_type=['N','A','D']
    for type in [0,1,2]:

        area=np.count_nonzero(dataland[:,1:]==type,axis=1).transpose()
        axs[1].plot(dataland[:,0],area,color=colors[type], label=land_type[type])

    axs[0].set_ylabel("Population size")
    axs[1].set_ylabel("Land area")
    axs[1].set_xlabel("Time")

    axs[1].legend()
    axs[0].set_title("a="+str(a)+", m="+str(m))
    plt.savefig("time_dynamics.jpg", dpi=500)
    plt.show()


    ## figure 2 mean ecosystem service provision with std and fragment size
    ################################################################################
    # fig, axs = plt.subplots(2,1,sharex='col')
    # axs[0].plot(dataclus[:,0],dataclus[:,5],color='tab:green')
    # axs[0].fill_between(dataclus[:,0], dataclus[:,5]+dataclus[:,6], dataclus[:,5]-dataclus[:,6], facecolor='tab:green', alpha=0.4)
    #
    # fragment_area_diff = np.abs(dataclus[:-1,3]-dataclus[1:,3])
    # ix2=np.where(fragment_area_diff==np.max(fragment_area_diff))[0][0]
    #
    # axs[1].scatter(dataclus[:ix2-25:25,0],dataclus[:ix2-25:25,2],s=dataclus[:ix2-25:25,3]*0.05,alpha=0.5,linewidths=2,color='tab:blue')
    # axs[1].scatter(dataclus[ix2:,0],dataclus[ix2:,2],s=dataclus[ix2:,3]*0.05,alpha=0.5,linewidths=2,color='tab:blue',label='Max area')
    #
    # axs[0].set_ylabel("ES provision")
    # axs[1].set_ylabel("Fragments")
    # axs[1].set_xlabel("Time")
    #
    # axs[1].legend()
    # axs[0].set_title("Percolation transition")
    # plt.savefig("perco_es.jpg", dpi=500)
    # plt.show()

    ## figure 3 patch size distribution before, at and after criticality
    ################################################################################

    # first identify where the critical transition happens, in order to do
    # that identify the largest jump in mean_fragment_area
    #
    # fig, axs = plt.subplots(1,3,sharex='row')
    #
    # # i need another datafile with the proper cluster information to acces the distribution
    # fragments_info=np.loadtxt("histpre.txt")
    # axs[0].hist(fragments_info[1:], 20, density=True, facecolor='tab:green', alpha=0.75)
    # plt.xscale('log')
    # fragments_info=np.loadtxt("histcrit.txt")
    # axs[1].hist(fragments_info[1:], 20, density=True, facecolor='tab:green', alpha=0.75)
    # plt.xscale('log')
    # fragments_info=np.loadtxt("histpost.txt")
    # axs[2].hist(fragments_info[1:], 20, density=True, facecolor='tab:green', alpha=0.75)
    # plt.xscale('log')
    # plt.show()
    # # fig.xscale("log");
    # # fig.xscale("log");
    # # fig.xscale("log");


plt.close()
