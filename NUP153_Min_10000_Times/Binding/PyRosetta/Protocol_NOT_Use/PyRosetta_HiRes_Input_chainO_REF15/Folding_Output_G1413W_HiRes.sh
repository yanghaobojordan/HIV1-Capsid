#!/bin/bash
#SBATCH -n 2
#SBATCH -t 01-00:00:00
#SBATCH --account=brubenst-condo

python Folding_Output_G1413W.py Folding_Output_G1413W.pdb
