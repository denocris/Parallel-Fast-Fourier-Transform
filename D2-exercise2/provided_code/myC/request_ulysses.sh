#!/bin/bash
#PBS -l nodes=16:ppn=20 -q wide -l walltime=00:30:00
#PBS -N fft

module load openmpi/1.8.3/gnu/4.9.2
module load fftw/3.3.4/gnu/4.9.2

cd /home/cdenobi/P2.8_seed/D2-exercise2/provided_code/myC
#cd $PBS_O_WORKDIR

make clean
make

num_nodes1="1 2 4 8"
#sizes="128"
#outfile="log_strong$size.txt"

#for size in $sizes; do
for node in $num_nodes1; do
    nproc=$(($node * 16))
    echo "------------- NUMBER OF PROCESSOR = $nproc -----------------" >> log_strong128.txt
    mpirun -np $nproc --map-by ppr:16:node:pe=1 ./diffusion.x 128 >> log_strong128.txt
  done
#done

num_nodes2="1 2 4 8 16"
#sizes="128"
#outfile="log_strong$size.txt"

#for size in $sizes; do
for node in $num_nodes2; do
    nproc=$(($node * 16))
    echo "------------- NUMBER OF PROCESSOR = $nproc -----------------" >> log_strong512.txt
    mpirun -np $nproc --map-by ppr:16:node:pe=1 ./diffusion.x 512 >> log_strong512.txt
  done
#done


num_nodes3="16"

#outfile="log_strong$size.txt"


#for size in $sizes; do
for node in $num_nodes3; do
    nproc=$(($node * 16))
    echo "------------- NUMBER OF PROCESSOR = $nproc -----------------" >> log_strong2048.txt
    mpirun -np $nproc --map-by ppr:16:node:pe=1 ./diffusion.x 2048 >> log_strong2048.txt
done
#done
