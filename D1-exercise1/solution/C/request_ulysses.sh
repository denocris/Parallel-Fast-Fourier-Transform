#!/bin/bash
#PBS -l nodes=16:ppn=20 -q wide -l walltime=00:10:00
#PBS -N fft

module load openmpi
module load fftw/3.3.4

cd /home/cdenobi/P2.8_seed/D1-exercise1/solution/C
#cd $PBS_O_WORKDIR

make

num_nodes="1 2 4 8 16"
dimension="128"


outfile="log_strong$size.txt"


for size in $dimension; do
  for node in $num_nodes; do
    nproc=$(($node * 8))
    echo "------------- NUMBER OF PROCESSOR = $nproc -----------------" >> $outfile
    mpirun -np $nproc --map-by ppr:8:node:pe=1 ./diffusion.x $size >> $outfile
  done
done
