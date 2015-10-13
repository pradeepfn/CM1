Pre-requisites

sudo mount -t tmpfs -o size=4G tmpfs /mnt/ramdisk
sudo chown -R username /mnt/ramdisk


##Running CM1###

set the 

nodex 
nodey

params in namelist.input.orig  file.

mpirun -np (nodex * nodey) ./cm1



Important:

the program reads the input file from every node.
