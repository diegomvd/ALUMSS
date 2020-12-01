#!/bin/bash

<<COMMENT1
	Script to run several simulations.

	Parameters are:

			1   2   3   4    5   6   7   8     9   10   11  12  13  14   15   16  17  18    19   20  21     22
			T  dtp  n  ao0  ai0  c0  r0  ys0  yn0  ess  a   w   g   Tao  Tai  Tr  Td  Tres  rho  d  dtsave  exp

COMMENT1

#run



./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 1 100 50 2 50 1 0 1 1 norest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 1 100 50 2 50 1 0 1 1 norest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 1 100 50 2 50 1 0 1 1 norest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 1 100 50 2 50 1 0 1 1 norest2&

./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 0.5 0 1 100 50 2 50 1 0 1 1 norest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 0.5 0 1 100 50 2 50 1 0 1 1 norest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 0.5 0 1 100 50 2 50 1 0 1 1 norest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 0.5 0 1 100 50 2 50 1 0 1 1 norest2&

./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 0.1 100 50 2 50 1 0 1 1 norest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 0.1 100 50 2 50 1 0 1 1 norest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 0.1 100 50 2 50 1 0 1 1 norest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 0.1 100 50 2 50 1 0 1 1 norest2&

./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 0.5 0 0.1 100 50 2 50 1 0 1 1 norest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 0.5 0 0.1 100 50 2 50 1 0 1 1 norest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 0.5 0 0.1 100 50 2 50 1 0 1 1 norest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 0.5 0 0.1 100 50 2 50 1 0 1 1 norest2&


./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 1 100 50 2 50 100 0.25 1 1 rest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 1 100 50 2 50 100 0.5 1 1 rest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 1 100 50 2 50 100 1.0 1 1 rest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 1 100 50 2 50 100 2.0 1 1 rest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 1 100 50 2 50 100 4.0 1 1 rest2&

./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 1 100 50 2 50 50 0.25 1 1 rest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 1 100 50 2 50 50 0.5 1 1 rest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 1 100 50 2 50 50 1.0 1 1 rest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 1 100 50 2 50 50 2.0 1 1 rest2&
./gillespie-ses 3000 0.1 40 0.1 0.0 1 1 10 10 1 1.0 0 1 100 50 2 50 50 4.0 1 1 rest2&
