#!/bin/bash
#SBATCH -n 2
#SBATCH -t 01-00:00:00
#SBATCH --account=brubenst-condo

python T1416M.py Folding_WT_8.pdb
