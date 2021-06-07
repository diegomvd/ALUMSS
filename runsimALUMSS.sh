#!/bin/bash

<<COMMENT1
	Script to run several simulations by wirting the execution commands (see below the parameter description)

	Parameters are:

			              1   2   3   4    5    6  7    8    9      10     11        12       13        14     15    16     17
names in code:	  	T  dtp  n  a0   d0  ksi  y0  sar   a      w     Tag        Tab      Tr        Td     d   dtsave  seed
names in paper:													y1	      z  alpha  omega  1/sigma   1/rho_L  1/rho_R	 1/rho_D

T = total simulation time
dtp = timestep to solve the population ODE
n = side length of the landscape in number of cells, total number of cells is n*n
a0 = initial fraction of agricultural land
d0 = initial fraction of degraded land
ksi = productivity of intensive agriculture
y0 = baseline productivity of low-intesive agriculture
sar = saturation exponent of ecosystem services with area of natural fragment
a = preference for agricultural intensification
w = clustering parameter of agricultural cells
Tag = decay rate of the average time until next expansion or intensification transition with respect to resource deficit (inverse of sensitivity to resource deficit)
Tab = growth rate of the average time until feritlity loss with respect to ecosystem service provision (inverse of fertility loss sensitivity to ecosysetm service provision)
Tr = decay rate of the average time until land recovery with respect to ecosystem service provision (inverse of recovery sensitivity to ecosystem service provision)
Td = growth rate of the average time until land degradation with respect to ecosystem service provision (inverse of degradation sensitivity to ecosysetm service provision)
d = distance in cell units over which ecosystem services are perceived, by default=1 in this study
dtsave = time step for saving data
seed = seed to initialize the random number generator

COMMENT1

# execution command

./alumss-exec 1000 0.1 40 0.2 0.0 1.2 0.2 0.25 0.0 0.0 0.1 50.0 5.0 50.0 1 1 186611281&
