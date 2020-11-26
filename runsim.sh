#!/bin/bash

<<COMMENT1
	Script to run several simulations.

	Parameters are:

			1   2   3   4    5   6   7   8     9   10   11  12  13  14   15   16  17  18    19   20  21     22
			T  dtp  n  ao0  ai0  c0  r0  ys0  yn0  ess  a   w   g   Tao  Tai  Tr  Td  Tres  rho  d  dtsave  exp

COMMENT1

#run



./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 1000 0 1 1 test0&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 0.5 0 10 80 20 2 50 1000 0 1 1 test0&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 0.0 0 10 80 20 2 50 1000 0 1 1 test0&

./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 40 2 50 1000 0 1 1 test1&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 0.5 0 10 80 40 2 50 1000 0 1 1 test1&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 0.0 0 10 80 40 2 50 1000 0 1 1 test1&

./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 100 50 2 50 1000 0 1 1 test2&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 0.5 0 10 100 50 2 50 1000 0 1 1 test2&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 0.0 0 10 100 50 2 50 1000 0 1 1 test2&

./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 1000 0 2 1 testd&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 1000 0 3 1 testd&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 1000 0 4 1 testd&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 1000 0 5 1 testd&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 1000 0 6 1 testd&

./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 25 1 1 1 testres&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 50 1 1 1 testres&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 100 1 1 1 testres&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 200 1 1 1 testres&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 400 1 1 1 testres&

./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 100 0.25 1 1 testres&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 100 0.5 1 1 testres&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 100 1 1 1 testres&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 100 2 1 1 testres&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 1.0 0 10 80 20 2 50 100 4 1 1 testres&

./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 0.5 0 10 80 20 2 50 100 0.25 1 1 testres&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 0.5 0 10 80 20 2 50 100 0.5 1 1 testres&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 0.5 0 10 80 20 2 50 100 1 1 1 testres&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 0.5 0 10 80 20 2 50 100 2 1 1 testres&
./gillespie-ses 3000 0.1 40 0.2 0.1 1 1 10 10 1 0.5 0 10 80 20 2 50 100 4 1 1 testres&
