#!/bin/bash
#BSUB -J comsol_far_does_in_fact_not_run_comsol_its_just_in_the_name_no_worries # name
#BSUB -o outfiles/far_field_comsol_test_%J.out # output file
#BSUB -q gpuh100
#BSUB -n 4 ## cores
#BSUB -R "rusage[mem=1GB]" 
#BSUB -W 4:00 # useable time in minutes
##BSUB -N # send mail when done
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=1:mode=exclusive_process"

module load gcc/12.3.0-binutils-2.40
module load cuda/12.2.2
# "Retinin2x2" or "demoleus2x2"
LAMBDA="325 425 525 625"
ANGLE="0 30 60 90"
for lambda in $LAMBDA;
do  
    for angle in $ANGLE;
    do
        python computeFarFieldPattern.py 2 $angle demoleus2x2 $lambda

        python computeFarFieldPattern.py 2 $angle Retinin2x2 $lambda
        
    done
done


exit 0