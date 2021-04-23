#!/bin/bash
#SBATCH -n 2
#SBATCH -t 02-00:00:00
#SBATCH --account=brubenst-condo

python WT_7.py hex-nup.pdb
