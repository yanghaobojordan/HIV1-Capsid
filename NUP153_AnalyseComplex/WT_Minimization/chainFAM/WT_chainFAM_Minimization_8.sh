#!/bin/bash
#SBATCH -n 2
#SBATCH -t 01-00:00:00
#SBATCH --account=brubenst-condo

python WT_chainFAM_Minimization_8.py hex-nup_chainFAM.pdb
