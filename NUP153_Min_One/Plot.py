# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 19:28:49 2020
@author: hyang
"""
import math
import numpy as np
from scipy.stats import reciprocal
import matplotlib.pyplot as plt
import pylab
from itertools import groupby
from statistics import mean
from numpy import array
from scipy import stats
import seaborn as sns
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
from scipy import optimize
sns.set(style="darkgrid")

def main():
    N=['P1411Y','P1411M','S1412P','S1412M','G1413W','G1413M','V1414W','V1414I','F1415G','F1415M','T1416R','T1416M','F1417G','F1417Y','G1418Y','G1418A']
    abs_pull=[12.92623273,32.13517882,63.03431965,50.11425081,57.41375282,86.71913347,20.94936705,16.68319625,28.36939488,78.7695472,7.80929512,35.41003671,26.50416747,65.38334736,17.91404869,42.47824555]
    norm=[52.0309552,67.88333333,66.10666667,64.06,69.22333333,60.7275,88.17666667,63.87333333,79.22333333,52.00666667,86.15666667,86.8,85.52666667,77.86333333,97.84333333,78.97333333]
    Fold=[206.5373612,86.4319566,79.01203002,10.65869171,116.7902131,69.8093994,337.9718907,56.53006166,31.86558884,14.46449091,3.23345021,-4.82628589,82.87364722,27.66255546,159.3639997,90.86094437]    
    
    fig = plt.figure(num=None, figsize=(11, 9),facecolor='w', edgecolor='k')   
    ax = fig.add_subplot(111)
    ax.set_xlabel("$\Delta\Delta$$G_{fold}$(kcal/mol)(PyRosetta)", size=15)
    ax.set_ylabel("TRIM-NUP153 Pulled Down By Hexamer Tubes (%)", size=15)
    #ax.set_ylabel("GFP+ Cells Resulting From TRIM Restriction(%)", size=15)
    ax.scatter(Fold, abs_pull, c='blue', marker='s', s=25)
    
    for x, y, n in zip(Fold, abs_pull, N): 
        plt.annotate(n, xy=(x,y), xytext=(17, -27),
                     textcoords='offset points',
                     ha="right", va='bottom',
                     bbox=dict(boxstyle='round,pad=0.5',
                               fc='yellow', alpha=0.5),
                               arrowprops=dict(arrowstyle = '->',connectionstyle='arc3,rad=0',color='black'))
               
main()