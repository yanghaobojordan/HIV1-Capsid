#!/bin/bash
#SBATCH -n 2
#SBATCH -t 01-00:00:00
#SBATCH --account=brubenst-condo

python Folding_Output_P1411Y.py Folding_Output_P1411Y.pdb
