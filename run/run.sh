#!/bin/bash


node_x=2
node_y=2

let no_of_procs=$node_x*$node_y
restart_freq=10

echo $node_x
echo $node_y

make rclean
cp namelist.input.orig namelist.input 
sed -i "s/nodex        =       1/nodex =  $node_x/" namelist.input
sed -i "s/nodey        =       1/nodey =  $node_y/" namelist.input
sed -i "s/rstfrq = -3600.0/rstfrq = $restart_freq/" namelist.input
mpirun -np $no_of_procs ./cm1.exe
