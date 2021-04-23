#!/bin/bash
#SBATCH -n 2
#SBATCH -t 01-00:00:00
#SBATCH --account=brubenst-condo

python S1412P_chainBCP.py WT_chainBCP_Minimization_2.pdb
