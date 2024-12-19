#!/bin/bash
#BSUB -J grid_far # name
#BSUB -o outfiles/grid_far_%J.out # output file
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
GRID="50 100 150 200 250 300 350 400 450 500 550 600 650 700 1000"
for grid in $GRID;
do
    for lambda in $LAMBDA;
    do  
        for angle in $ANGLE;
        do
            python computeFarFieldPattern.py 4 $angle demoleus2x2 $lambda $grid

            python computeFarFieldPattern.py 4 $angle Retinin2x2 $lambda $grid
            
        done
    done
done

exit 0