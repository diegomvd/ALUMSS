#!/bin/bash

<<COMMENT1
	Script to run several simulations.

	Parameters are:

      1    2   3  4    5    6  7   8  9  10  11  12   13   14  15   16   
      T   dtp  n  p0  cg0  r0  y0  a  w  m   g   Ta   phi  Tr  Td  dtsave

COMMENT1

./gillespie-ses 3000 0.1 40 5 1 1 5 0.0 0 1 0.1 1 3 1.5 50 1 &
./gillespie-ses 3000 0.1 40 5 1 1 5 0.25 0 1 0.1 1 3 1.5 50 1 &
./gillespie-ses 3000 0.1 40 5 1 1 5 0.5 0 1 0.1 1 3 1.5 50 1 &
./gillespie-ses 3000 0.1 40 5 1 1 5 1.0 0 1 0.1 1 3 1.5 50 1 &
./gillespie-ses 3000 0.1 40 5 1 1 5 1.5 0 1 0.1 1 3 1.5 50 1 &
./gillespie-ses 3000 0.1 40 5 1 1 5 2.0 0 1 0.1 1 3 1.5 50 1 &
./gillespie-ses 3000 0.1 40 5 1 1 5 4.0 0 1 0.1 1 3 1.5 50 1 &
