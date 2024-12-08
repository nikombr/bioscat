#!/bin/bash
#BSUB -J segment_test # name
#BSUB -o outfiles/segment_test_%J.out # output file
#BSUB -q gpuh100
#BSUB -n 64 ## cores
#BSUB -R "rusage[mem=1GB]" 
#BSUB -W 24:00 # useable time in minutes
##BSUB -N # send mail when done
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=1:mode=exclusive_process"

module load gcc/12.3.0-binutils-2.40
module load cuda/12.2.2
# "Retinin2x2" or "demoleus2x2"

OMP_NUM_THREADS=64 python forward.py 1 1 0 demoleus2x2

OMP_NUM_THREADS=64 python forward.py 1 2 0 demoleus2x2

OMP_NUM_THREADS=64 python forward.py 1 3 0 demoleus2x2

OMP_NUM_THREADS=64 python forward.py 1 4 0 demoleus2x2

OMP_NUM_THREADS=64 python forward.py 1 5 0 demoleus2x2

OMP_NUM_THREADS=64 python forward.py 1 1 90 demoleus2x2

OMP_NUM_THREADS=64 python forward.py 1 2 90 demoleus2x2

OMP_NUM_THREADS=64 python forward.py 1 3 90 demoleus2x2

OMP_NUM_THREADS=64 python forward.py 1 4 90 demoleus2x2

OMP_NUM_THREADS=64 python forward.py 1 5 90 demoleus2x2

OMP_NUM_THREADS=64 python forward.py 1 1 0 Retinin2x2

OMP_NUM_THREADS=64 python forward.py 1 2 0 Retinin2x2

OMP_NUM_THREADS=64 python forward.py 1 3 0 Retinin2x2

OMP_NUM_THREADS=64 python forward.py 1 4 0 Retinin2x2

OMP_NUM_THREADS=64 python forward.py 1 5 0 Retinin2x2

OMP_NUM_THREADS=64 python forward.py 1 1 90 Retinin2x2

OMP_NUM_THREADS=64 python forward.py 1 2 90 Retinin2x2

OMP_NUM_THREADS=64 python forward.py 1 3 90 Retinin2x2

OMP_NUM_THREADS=64 python forward.py 1 4 90 Retinin2x2

OMP_NUM_THREADS=64 python forward.py 1 5 90 Retinin2x2