#!/bin/bash
#SBATCH -n 2
#SBATCH -t 01-00:00:00
#SBATCH --account=brubenst-condo

python WT_chainEFV_Minimization_6.py hex-nup_chainEFV.pdb
