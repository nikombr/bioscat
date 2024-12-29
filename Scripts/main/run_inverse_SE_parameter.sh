#!/bin/bash
#BSUB -J inverseSE # name
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

num_execution=2
chainNum=1
echo $num_execution

if [ $num_execution == 1 ];
then
    python inverse.py 2 $chainNum Retinin2x2 angle_resolved 0.7 0.07
    python inverse.py 2 $chainNum Retinin2x2 angle_resolved 0.7 0.1
    python inverse.py 2 $chainNum Retinin2x2 angle_resolved 0.7 0.04
    python inverse.py 2 $chainNum Retinin2x2 angle_resolved 0.5 0.07
    python inverse.py 2 $chainNum Retinin2x2 angle_resolved 0.9 0.07
fi

if [ $num_execution == 2 ];
then
    python inverse.py 2 $chainNum demoleus2x2 angle_resolved 0.7 0.07
    python inverse.py 2 $chainNum demoleus2x2 angle_resolved 0.7 0.1
    python inverse.py 2 $chainNum demoleus2x2 angle_resolved 0.7 0.04
    python inverse.py 2 $chainNum demoleus2x2 angle_resolved 0.5 0.07
    python inverse.py 2 $chainNum demoleus2x2 angle_resolved 0.9 0.07
fi

exit

if [ $num_execution == 1 ];
then
    python inverse.py 2 $chainNum demoleus2x2 angle_resolved 0.7 0.07
    python inverse.py 2 $chainNum Retinin2x2 angle_resolved 0.7 0.07
fi

if [ $num_execution == 2 ];
then
    python inverse.py 2 $chainNum demoleus2x2 angle_resolved 0.7 0.1
    python inverse.py 2 $chainNum Retinin2x2 angle_resolved 0.7 0.1
fi

if [ $num_execution == 3 ];
then
    python inverse.py 2 $chainNum demoleus2x2 angle_resolved 0.7 0.04
    python inverse.py 2 $chainNum Retinin2x2 angle_resolved 0.7 0.04
fi

if [ $num_execution == 4 ];
then
    python inverse.py 2 $chainNum demoleus2x2 angle_resolved 0.5 0.07
    python inverse.py 2 $chainNum Retinin2x2 angle_resolved 0.5 0.07
fi
if [ $num_execution == 5 ];
then
    python inverse.py 2 $chainNum demoleus2x2 angle_resolved 0.9 0.07
    python inverse.py 2 $chainNum Retinin2x2 angle_resolved 0.9 0.07
fi

## ELL="0.04 0.07 0.1";
## TAU="0.5 0.7 0.9";

## for ell in $ELL;
## do
##     for tau in $TAU;
##     do
##         OMP_NUM_THREADS=64 python inverse.py 2 demoleus2x2 angle_resolved $tau $ell
##     done
## done