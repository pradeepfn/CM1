#!/bin/bash
echo 0 >| notify/cm1.notify

node_x=4
node_y=1

let no_of_procs=$node_x*$node_y
restart_freq=10

echo $node_x
echo $node_y

make restartclean
cp namelist.input.orig namelist.input 
#cp phoenix.config.run phoenix.config
sed -i "s/nodex        =       1/nodex =  $node_x/" namelist.input
sed -i "s/nodey        =       1/nodey =  $node_y/" namelist.input
sed -i "s/rstfrq = -3600.0/rstfrq = $restart_freq/" namelist.input
#mpirun -np $no_of_procs -f host_file ./cm1.exe
mpirun -np $no_of_procs --bind-to core ../../phoenix/bin/mpiformat /dev/shm 1000
mpirun -np $no_of_procs --bind-to core ./cm1.exe
