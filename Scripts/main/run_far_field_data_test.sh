#!/bin/bash
#BSUB -J datafar # name
#BSUB -o outfiles/far_field_data_test_%J.out # output file
#BSUB -q gpuh100
#BSUB -n 4 ## cores
#BSUB -R "rusage[mem=1GB]" 
#BSUB -W 24:00 # useable time in minutes
##BSUB -N # send mail when done
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=1:mode=exclusive_process"

module load gcc/12.3.0-binutils-2.40
module load cuda/12.2.2
# "Retinin2x2" or "demoleus2x2"

python computeFarFieldPattern.py 5 flattest
## python computeFarFieldPattern.py 5 demoleus2x2

## python computeFarFieldPattern.py 5 Retinin2x2
    



exit 0