#!/bin/bash
#SBATCH -n 2
#SBATCH -t 02-00:00:00
#SBATCH --account=brubenst-condo

python HiRes_Cross_Dock.py Folding_WT_chainABP_8.pdb
