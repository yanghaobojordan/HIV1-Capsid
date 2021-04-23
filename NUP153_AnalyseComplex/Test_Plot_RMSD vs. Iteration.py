# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 12:57:23 2020

@author: hyang
"""

import seaborn as sns
import matplotlib.patches as mpatches
import collections
import random
import numpy as np
import matplotlib.pyplot as plt
import pylab
from statistics import mean
from numpy import array
from scipy import stats
from matplotlib.patches import Rectangle
sns.set(style="darkgrid")

def main():

    file=open("Test_Plot_RMSD vs. Iteration.txt","r")
    Iteration=[]
    RMSD=[]
    Energy=[]
    for i in file:
        i=i.rstrip().split()
        Iteration.append(float(i[1]))
        Energy.append(float(i[0]))
        RMSD.append(float(i[2]))
    file.close()
    fig = plt.figure(num=None, figsize=(10, 8),facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    ax.scatter(Iteration, Energy, c='dodgerblue', marker='o', s=4.0)
    ax.set_xlabel("Iteration", size=17)
    ax.set_ylabel("RMSD", size=17)
    #plt.ylim(0,1)
main()