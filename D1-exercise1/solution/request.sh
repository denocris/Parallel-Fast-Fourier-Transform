#!/bin/bash

#PBS -N ExFftw
#PBS -l walltime=00:20:30
#PBS -l nodes=4:ppn=20

module load openmpi/1.8.3/gnu/4.9.2
module load fftw/3.3.4/gnu/4.9.2

for nprocs in 1 2 4 8 16 20;
do
  /usr/bin/time -p mpirun -np $nprocs ./diffusion.x in 2>&1 timing.dat
done

cat timing.dat | grep real | awk '{print $2}' > timing_cleaned.dat
