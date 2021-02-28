#!/bin/bash
#BSUB -J Merfish_450FOV
#BSUB -n 12
#BSUB -R span[ptile=2]
#BSUB -R rusage[mem=4]
#BSUB -W 0:20

cd $LS_SUBCWD

python start_file.py
