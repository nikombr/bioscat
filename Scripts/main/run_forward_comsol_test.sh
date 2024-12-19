#!/bin/bash
#BSUB -J comsol_forward_does_in_fact_not_run_comsol_its_just_in_the_name_no_worries # name
#BSUB -o outfiles/comsol_test_forward_%J.out # output file
#BSUB -q gpuh100
#BSUB -n 4 ## cores
#BSUB -R "rusage[mem=8GB]" 
#BSUB -W 4:00 # useable time in minutes
##BSUB -N # send mail when done
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=1:mode=exclusive_process"

module load gcc/12.3.0-binutils-2.40
module load cuda/12.2.2
# "Retinin2x2" or "demoleus2x2"
ANGLE="0 30 60 90"

for angle in $ANGLE;
do
echo $angle
python forward.py 2 $angle demoleus2x2 300

python forward.py 2 $angle Retinin2x2 300

done
