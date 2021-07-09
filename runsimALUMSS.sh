#!/bin/bash

<<COMMENT1
	Script to run several simulations by wirting the execution commands (see below the parameter description)

	Parameters are:

			              1   2   3   4    5    6   7    8    9   10    11   12    13    14   15   16     17
names in code:	  	T   n  a0  d0   nF    a  mS   wS    z   dES   y0   y1    sFL   sR   sD  dtSave  seed

T = total simulation time
n = side length of the landscape in number of cells, total number of cells is n*n
a0 = initial fraction of agricultural land
d0 = initial fraction of degraded land
nF = number of farms
a = fraction of sparing farms
mS = total sensitivity to resource demand
wS = width of the sensitivity spectrum from 0 to 1 : std(sA_i) ~ wS
z = ES-Area relationship saturation exponent
dES = ecosystemic characteristic spatial scale (connection distance and ES flow distance)
y0 = low intense baseline productivity
y1 = high intense productivity
sFL = fertility loss sensitivity to ES provision
sR = land recovery sensitivity to ES provision
sD = land degradation sensitivity to ES provision
dtSave = time step for saving data
seed = seed to initialize the random number generator

COMMENT1

# execution command

./alumss-exec 10 10 0.2 0.0 1.0 0.0 10 0.0 0.25 1.1 0.2 1.2 0.02 0.2 0.02 1 186611281&
