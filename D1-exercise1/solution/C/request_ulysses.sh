#!/bin/bash
#PBS -l nodes=16:ppn=20 -q wide -l walltime=08:00:00

module load openmpi
module load fftw/3.3.4

cd P2.8_seed/D1-exercise1/solution/C
#cd $PBS_O_WORKDIR

num_nodes="1 2 4 8 16"
dimension="64"

for size in $dimension; do
  for node in $num_nodes; do
    nproc=$(($node * 16))
    mpirun -np $nproc --map-by ppr:16:node:pe=1 ./diffusion.x $size >> strong_$size.txt
  done
done
