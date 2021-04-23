#!/bin/bash
#SBATCH -n 2
#SBATCH -t 01-00:00:00
#SBATCH --account=brubenst-condo

python WT_chainCDT_Minimization_5.py hex-nup_chainCDT.pdb
