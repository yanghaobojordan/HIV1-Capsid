#!/bin/bash
#SBATCH -n 2
#SBATCH -t 01-00:00:00
#SBATCH --account=brubenst-condo

python WT_chainCDT_Minimization_1.py hex-nup_chainCDT.pdb
