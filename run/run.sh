#!/bin/bash

make restartclean
cp namelist.input.orig namelist.input 
sed -i "s/rstfrq = -3600.0/rstfrq = 10/" namelist.input
./cm1.exe
