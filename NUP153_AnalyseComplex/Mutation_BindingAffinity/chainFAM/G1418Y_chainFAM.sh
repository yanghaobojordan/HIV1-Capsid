#!/bin/bash
#SBATCH -n 2
#SBATCH -t 01-00:00:00
#SBATCH --account=brubenst-condo

python G1418Y_chainFAM.py WT_chainFAM_Minimization_7.pdb
