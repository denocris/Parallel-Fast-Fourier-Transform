#!/bin/bash

#PBS -N fft
#PBS -l walltime=00:30:00
#PBS -l nodes=4:ppn=20

module load openmpi
module load fftw/3.3.4

cd /home/cdenobi/P2.8_seed/D1-exercise1/solution/C

#cd $PBS_O_WORKDIR

for NPE in 80 #1 2 4 8 16 20 #32 40
do
    for COUNT in 1 2 3 4
    do
	mpirun -np $NPE diffusion.x
    done
done
