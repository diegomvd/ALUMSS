#!/bin/bash

<<COMMENT1
	Script to run several simulations.

	Parameters are:

      1    2   3  4    5    6  7   8  9  10  11  12   13   14  15   16		17	18
      T   dtp  n  a0  cg0  r0  y0  a  w  m   g   Ta   phi  Tr  Td  dtsave  k  Ti

COMMENT1

# test local maintenance in new degradation model

./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.2 0 1.1 0.1 2 3 2 50 1 10 10 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.3 0 1.1 0.1 2 3 2 50 1 10 10 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.4 0 1.1 0.1 2 3 2 50 1 10 10 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.5 0 1.1 0.1 2 3 2 50 1 10 10 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.6 0 1.1 0.1 2 3 2 50 1 10 10 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.7 0 1.1 0.1 2 3 2 50 1 10 10 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.8 0 1.1 0.1 2 3 2 50 1 10 10 &

./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.2 0 1.1 0.1 2 3 2 50 1 10 20 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.3 0 1.1 0.1 2 3 2 50 1 10 20 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.4 0 1.1 0.1 2 3 2 50 1 10 20 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.5 0 1.1 0.1 2 3 2 50 1 10 20 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.6 0 1.1 0.1 2 3 2 50 1 10 20 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.7 0 1.1 0.1 2 3 2 50 1 10 20 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.8 0 1.1 0.1 2 3 2 50 1 10 20 &

./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.2 0 1.1 0.1 2 3 2 50 1 10 40 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.3 0 1.1 0.1 2 3 2 50 1 10 40 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.4 0 1.1 0.1 2 3 2 50 1 10 40 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.5 0 1.1 0.1 2 3 2 50 1 10 40 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.6 0 1.1 0.1 2 3 2 50 1 10 40 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.7 0 1.1 0.1 2 3 2 50 1 10 40 &
./gillespie-ses 2000 0.1 40 0.10 1 1 5 0.8 0 1.1 0.1 2 3 2 50 1 10 40 &
