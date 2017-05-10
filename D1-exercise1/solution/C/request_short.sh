#!/bin/bash
#PBS -l nodes=1:ppn=20 -q regular -l walltime=00:01:00
module load openmpi/1.8.3/gnu/4.9.2
module load fftw/3.3.4/gnu/4.9.2

cd /home/cdenobi/P2.8_seed/D1-exercise1/solution/C

make

processes="1 2 4 8 16"
sizes="48"

for size in $sizes; do
outfile="log_strong$size.txt"
#echo "---------------------STRONG SCALABILITY---------------------" > $outfile
	for p in $processes; do
		#p=$(($node * 16))
		echo "-------------NUMBER OF PROCESSOR=$p-----------------" >> $outfile
		mpirun -np $p --map-by ppr:16:node:pe=1 ./diffusion.x $size >> $outfile

	done

done
