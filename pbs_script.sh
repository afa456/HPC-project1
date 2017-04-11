#!/bin/sh

#PBS -q class
#PBS -l nodes=4:sixcore
#PBS -l walltime=00:10:00
#PBS -N hpc-spring17-prog1-master

#Commands:
# change to our project directory
#cd $HOME/hpc-spring17-prog1-master
# hardcode MPI path
#MPIRUN=/usr/lib64/openmpi/bin/mpirun
# loop over number of processors (our 4 nodes job can run up to 48)
p=3
#do 
MPIRUN -np p ./poly-eval sample-constants.txt sample-values.txt
#done