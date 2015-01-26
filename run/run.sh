#!/bin/bash

make restartclean
cp namelist.input.orig namelist.input 
sed -i "s/rstfrq = -3600.0/rstfrq = 10/" namelist.input
mpirun -np 4 ./cm1.exe
