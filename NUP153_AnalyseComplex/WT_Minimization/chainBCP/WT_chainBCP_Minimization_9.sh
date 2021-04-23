#!/bin/bash
#SBATCH -n 2
#SBATCH -t 01-00:00:00
#SBATCH --account=brubenst-condo

python WT_chainBCP_Minimization_9.py hex-nup_chainBCP.pdb
