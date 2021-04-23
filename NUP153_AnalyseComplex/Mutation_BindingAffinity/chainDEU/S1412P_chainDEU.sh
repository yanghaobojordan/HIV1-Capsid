#!/bin/bash
#SBATCH -n 2
#SBATCH -t 01-00:00:00
#SBATCH --account=brubenst-condo

python S1412P_chainDEU.py WT_chainDEU_Minimization_2.pdb
