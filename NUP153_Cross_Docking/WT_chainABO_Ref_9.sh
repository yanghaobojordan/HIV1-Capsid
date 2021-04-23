#!/bin/bash
#SBATCH -n 2
#SBATCH -t 01-00:00:00
#SBATCH --account=brubenst-condo

python WT_chainABO_Ref_9.py hex-nup_ABO_Ref.pdb
