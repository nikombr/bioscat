#!/bin/bash
#BSUB -J segment_far # name
#BSUB -o outfiles/far_field_segment_test_%J.out # output file
#BSUB -q gpuh100
#BSUB -n 5 ## cores
#BSUB -R "rusage[mem=1GB]" 
#BSUB -W 4:00 # useable time in minutes
##BSUB -N # send mail when done
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=1:mode=exclusive_process"

module load gcc/12.3.0-binutils-2.40
module load cuda/12.2.2
# "Retinin2x2" or "demoleus2x2"

ANGLE="0 90"
N_SEG="1 2 3 4 5"

for angle in $ANGLE;
do
    for n in $N_SEG;
    do
        OMP_NUM_THREADS=5 python computeFarFieldPattern.py 1 $n $angle demoleus2x2

        OMP_NUM_THREADS=5 python computeFarFieldPattern.py 1 $n $angle Retinin2x2
done

done
