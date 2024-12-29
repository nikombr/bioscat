#!/bin/bash
#BSUB -J prof->h01 # name
#BSUB -o outfiles/profile_forward_%J.out # output file
#BSUB -q gpuh100
#BSUB -n 8 ## cores
#BSUB -R "rusage[mem=8GB]" 
#BSUB -W 24:00 # useable time in minutes
##BSUB -N # send mail when done
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=1:mode=exclusive_process"

module load gcc/12.3.0-binutils-2.40
module load cuda/12.2.2

OMP_NUM_THREADS=1 python profileForward.py 2 1
