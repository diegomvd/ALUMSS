#!/bin/bash

<<COMMENT1
	Script to run several simulations.

	Parameters are:

			1   2   3   4    5    6    7   8   9   10    11    12  13  14   15     16
			T  dtp  n  ao0  ai0  ksi  sar  a   w   Tag   Tab   Tr  Td   d  dtsave  exp

COMMENT1

#run

./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 0.0 0 0.1 50 2 50 1 1 testa&
./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 0.1 0 0.1 50 2 50 1 1 testa&
./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 0.2 0 0.1 50 2 50 1 1 testa&
./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 0.4 0 0.1 50 2 50 1 1 testa&
./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 0.8 0 0.1 50 2 50 1 1 testa&
./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 1.6 0 0.1 50 2 50 1 1 testa&
./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 3.2 0 0.1 50 2 50 1 1 testa&
./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 6.4 0 0.1 50 2 50 1 1 testa&

./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 0.0 0 0.05 50 2 50 1 1 testb&
./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 0.1 0 0.05 50 2 50 1 1 testb&
./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 0.2 0 0.05 50 2 50 1 1 testb&
./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 0.4 0 0.05 50 2 50 1 1 testb&
./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 0.8 0 0.05 50 2 50 1 1 testb&
./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 1.6 0 0.05 50 2 50 1 1 testb&
./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 3.2 0 0.05 50 2 50 1 1 testb&
./gillespie-ses 4000 0.1 40 0.2 0 1 0.3 6.4 0 0.05 50 2 50 1 1 testb&
