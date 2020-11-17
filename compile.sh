#!/bin/bash

c++ gill_main.cpp gill_functions.cpp cokus3.c -o gillespie-ses-local -lgsl -lgslcblas -lm -Wall -Weffc++ --std=c++17 -lstdc++fs
