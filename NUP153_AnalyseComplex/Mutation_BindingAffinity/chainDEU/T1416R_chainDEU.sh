#!/bin/bash
#SBATCH -n 2
#SBATCH -t 01-00:00:00
#SBATCH --account=brubenst-condo

python T1416R_chainDEU.py WT_chainDEU_Minimization_2.pdb
