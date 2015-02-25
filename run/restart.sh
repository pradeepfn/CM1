#!/bin/bash
echo 1>| notify/cm1.notify
node_x=4
node_y=4

let no_of_procs=$node_x*$node_y
restart_freq=1000

echo $node_x
echo $node_y
cp namelist.input.orig namelist.input 
sed -i "s/nodex        =       1/nodex =  $node_x/" namelist.input
sed -i "s/nodey        =       1/nodey =  $node_y/" namelist.input
sed -i "s/rstfrq = -3600.0/rstfrq = $restart_freq/" namelist.input
sed -i "s/irst      =  0/irst      =  5/" namelist.input
sed -i "s/rstnum    =  1/irst      =  3/" namelist.input
mpirun -np $no_of_procs -f host_file ./cm1.exe
