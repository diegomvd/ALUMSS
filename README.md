Agricultural land-use management spatial simulation (ALUMSS)

This repository contains the code used to simulate the coupled evolution of a human population densitiy and the spatially explicit landscape the population manages and modifies in order to obtain agricultural resources. 

The model is described in depth in the article 

The code to run the simulations is written in C++ and is called mainALUMSS.cpp. The functions required to run the code are declared in functionsALUMSS.h and implemented in functionsALUMSS.cpp. Compilation can be done by executing the compileALUMSS.sh bash script to create the executable alumss-exec. An easy way to perform simulations is to execute runsimALUMSS.sh bash script where a list of simulations can be specified. A guide to specify the paramter values is found in the script.     

For more elaborated numerical experients and model explorations we used the OpenMOLE tool (https://openmole.org/). In the file explorationALUMSS.oms the simulation model is embedded in OpenMOLE and several numerical experiments are defined. To a guide on how to create new experiments refer to the OpenMOLE documentation (https://openmole.org/Documentation.html). The file sensitivityALUMSS.oms contains the implementation of a Morris Sensitivity Analysis made with OpenMOLE.

The programs to analyze and plot the outputs from the simulations are written in python. Here we provide all the programs used to plot the figures appearing in the main text and the supplementary materials of article. We also provide the figures from the main text in PDF format. The figure were produced with plotFig1-intensification.py, plotFig2-fragmentation.py, plotFig3-clustering.py and plotFig4-sharing-sparing.py. 

If you have any questions regarding the code from the simulations or the analysis or have trouble performing the simulations please feel free to email me at: 

diego.bengochea-paz@sete.cnrs.fr
or
bengo.bengo@pm.me
