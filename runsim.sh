#!/bin/bash

<<COMMENT1
	Script to run several simulations.

	Parameters are:

      1    2   3  4    5    6  7   8  9  10  11  12   13   14  15   16
      T   dtp  n  a0  cg0  r0  y0  a  w  m   g   Ta   phi  Tr  Td  dtsave

COMMENT1

# test local maintenance in new degradation model
./gillespie-ses 2000 0.1 40 0.10 1 1 5 1 0 1.1 1 2 3 2 50 1 &
./gillespie-ses 2000 0.1 40 0.15 1 1 5 1 0 1.1 1 2 3 2 50 1 &
./gillespie-ses 2000 0.1 40 0.20 1 1 5 1 0 1.1 1 2 3 2 50 1 &
./gillespie-ses 2000 0.1 40 0.25 1 1 5 1 0 1.1 1 2 3 2 50 1 &
./gillespie-ses 2000 0.1 40 0.30 1 1 5 1 0 1.1 1 2 3 2 50 1 &
./gillespie-ses 2000 0.1 40 0.35 1 1 5 1 0 1.1 1 2 3 2 50 1 &
./gillespie-ses 2000 0.1 40 0.40 1 1 5 1 0 1.1 1 2 3 2 50 1 &
./gillespie-ses 2000 0.1 40 0.45 1 1 5 1 0 1.1 1 2 3 2 50 1 &
./gillespie-ses 2000 0.1 40 0.50 1 1 5 1 0 1.1 1 2 3 2 50 1 &
./gillespie-ses 2000 0.1 40 0.55 1 1 5 1 0 1.1 1 2 3 2 50 1 &
./gillespie-ses 2000 0.1 40 0.60 1 1 5 1 0 1.1 1 2 3 2 50 1 &
./gillespie-ses 2000 0.1 40 0.65 1 1 5 1 0 1.1 1 2 3 2 50 1 &
./gillespie-ses 2000 0.1 40 0.70 1 1 5 1 0 1.1 1 2 3 2 50 1 &
