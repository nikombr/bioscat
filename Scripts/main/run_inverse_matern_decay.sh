#!/bin/bash
#BSUB -J invDecayM # name
#BSUB -o outfiles/inverse_%J.out # output file
#BSUB -q gpuh100
#BSUB -n 4 ## cores
#BSUB -R "rusage[mem=1GB]" 
#BSUB -W 24:00 # useable time in minutes
##BSUB -N # send mail when done
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=1:mode=exclusive_process"

module load gcc/12.3.0-binutils-2.40
module load cuda/12.2.2

num_execution=4
chainNum=1
echo $num_execution

if [ $num_execution == 1 ];
then
    python inverse.py 3 $chainNum Retinin2x2 angle_resolved 0.1 1e-3
    python inverse.py 3 $chainNum demoleus2x2 angle_resolved 0.1 1e-3
fi

if [ $num_execution == 2 ];
then
    python inverse.py 3 $chainNum Retinin2x2 angle_resolved 0.3 1e-3
    python inverse.py 3 $chainNum demoleus2x2 angle_resolved 0.3 1e-3
fi
if [ $num_execution == 3 ];
then
    python inverse.py 3 $chainNum Retinin2x2 angle_resolved 0.1 1e-2
    python inverse.py 3 $chainNum demoleus2x2 angle_resolved 0.1 1e-2
fi

if [ $num_execution == 4 ];
then
    python inverse.py 3 $chainNum Retinin2x2 angle_resolved 0.3 1e-2
    python inverse.py 3 $chainNum demoleus2x2 angle_resolved 0.3 1e-2
fi

