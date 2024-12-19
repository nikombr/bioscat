#!/bin/bash
#BSUB -J comsol # name
#BSUB -o outfiles/comsol_%J.out # output file
#BSUB -q hpc
#BSUB -n 64 ## cores
#BSUB -R "rusage[mem=1GB]" 
#BSUB -W 24:00 # useable time in minutes
##BSUB -N # send mail when done
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=1:mode=exclusive_process"

comsolmatlab -tmpdir work3/s194146/comsoltmp -recoverydir work3/s194146/comsolrecovery 

sleep 60

matlab -batch generateDataForInverse