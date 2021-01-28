#!/bin/bash

<<COMMENT1
	Script to run several simulations.

	Parameters are:

			1   2   3   4    5    6    7   8   9   10    11    12  13  14   15     16
			T  dtp  n  ao0  ai0  ksi  sar  a   w   Tag   Tab   Tr  Td   d  dtsave  exp

COMMENT1

#run

./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 0.0 0 0.05 50 5 50 1 1 a0&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 0.0 1 0.05 50 5 50 1 1 a0&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 0.0 2 0.05 50 5 50 1 1 a0&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 0.0 3 0.05 50 5 50 1 1 a0&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 0.0 4 0.05 50 5 50 1 1 a0&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 0.0 5 0.05 50 5 50 1 1 a0&

./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 0.5 0 0.05 50 5 50 1 1 a05&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 0.5 1 0.05 50 5 50 1 1 a05&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 0.5 2 0.05 50 5 50 1 1 a05&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 0.5 3 0.05 50 5 50 1 1 a05&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 0.5 4 0.05 50 5 50 1 1 a05&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 0.5 5 0.05 50 5 50 1 1 a05&

./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 1 0 0.05 50 5 50 1 1 a1&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 1 1 0.05 50 5 50 1 1 a1&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 1 2 0.05 50 5 50 1 1 a1&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 1 3 0.05 50 5 50 1 1 a1&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 1 4 0.05 50 5 50 1 1 a1&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 1 5 0.05 50 5 50 1 1 a1&

./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 2 0 0.05 50 5 50 1 1 a2&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 2 1 0.05 50 5 50 1 1 a2&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 2 2 0.05 50 5 50 1 1 a2&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 2 3 0.05 50 5 50 1 1 a2&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 2 4 0.05 50 5 50 1 1 a2&
./gillespie-ses 2000 0.1 40 0.2 0 1 0.3 2 5 0.05 50 5 50 1 1 a2&
