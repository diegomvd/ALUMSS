#!/bin/bash

<<COMMENT1
	Script to run several simulations.

	Parameters are:

      1    2   3  4    5    6  7   8  9  10  11  12   13   14  15   16
      T   dtp  n  a0  cg0  r0  y0  a  w  m   g   Ta   phi  Tr  Td  dtsave

COMMENT1

# test for effect of global maintenance
./gillespie-ses 2000 0.1 40 0.1 1 1 5 0.0 0 1.1 0.1 1 3 1.5 50 1 &
./gillespie-ses 2000 0.1 40 0.1 1 1 5 0.0 0 1.6 0.1 1 3 1.5 50 1 &
./gillespie-ses 2000 0.1 40 0.1 1 1 5 0.0 0 2.0 0.1 1 3 1.5 50 1 &

# # threshold on sustainability of cropped patch at 0.75 exposure to nature
# ./gillespie-ses 3000 0.1 40 0.1 1 1 5 0.0 0 1.6 0.1 1 3 1.5 50 1 &
# ./gillespie-ses 3000 0.1 40 0.1 1 1 5 0.1 0 1.6 0.1 1 3 1.5 50 1 &
# ./gillespie-ses 3000 0.1 40 0.1 1 1 5 0.2 0 1.6 0.1 1 3 1.5 50 1 &

# ./gillespie-ses 3000 0.1 40 1 1 1 5 4.0 0 1 0.1 1 3 1.5 50 1 &
