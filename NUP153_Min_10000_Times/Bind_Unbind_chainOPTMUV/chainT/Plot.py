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
from matplotlib.patches import Rectangle
sns.set(style="darkgrid")

def main():
    N_Disrupt=['P1411Y','S1412P','G1413W','V1414W','F1415G','T1416R','F1417G','G1418Y']
    abs_pull_Disrupt=[12.92623273,63.03431965,57.41375282,20.94936705,28.36939488,7.80929512,26.50416747,17.91404869]
    norm_Disrupt=[68.03,66.10666667,69.22333333,88.17666667,79.22333333,86.15666667,85.52666667,97.84333333]
    FoldX_Disrupt=[3.2913,1.413089,0.5442,15.19362,2.597332,0.420218,5.577509,7.891621]
    PyRosetta_Disrupt=[8.43477873,5.36092693,2.92390854,23.56816268,5.65028464,0.39317595,11.5967968,1.15309509]
    
    N_WT_like=['P1411M','S1412M','G1413M','V1414I','F1415M','T1416M','F1417Y','G1418A']
    abs_pull_WT_like=[32.13517882,50.11425081,86.71913347,16.68319625,78.7695472,35.41003671,65.38334736,42.47824555]
    norm_WT_like=[67.88333333,64.06,60.7275,63.87333333,52.00666667,86.8,77.86333333,78.97333333]
    FoldX_WT_like=[-0.083443,-0.587499,-0.757159,0.260271,0.173853,-2.097057,0.937911,1.725118]
    PyRosetta_WT_like=[2.48188082,-1.14332058,3.12628322,2.12310464,3.57530799,-0.18906652,7.1868399,-0.13876031]

    fig = plt.figure(num=None, figsize=(11, 9),facecolor='w', edgecolor='k')   
    ax = fig.add_subplot(111)
    ax.set_xlabel("$\Delta\Delta$$G_{bind}$(kcal/mol)(PyRosetta)(chainO)", size=15)
    ax.set_ylabel("$\Delta\Delta$$G_{bind}$(kcal/mol)(FoldX)", size=15)
    #ax.set_ylabel("TRIM-NUP153 Pulled Down By Hexamer Tubes (%)", size=15)
    #ax.set_ylabel("GFP+ Cells Resulting From TRIM Restriction(%)", size=15)
    ax.scatter(PyRosetta_Disrupt, FoldX_Disrupt, c='red', marker='s', s=25)
    ax.scatter(PyRosetta_WT_like, FoldX_WT_like, c='blue', marker='s', s=25)
    #ax.errorbar(PyRosetta_Disrupt, norm_Disrupt, xerr=PyRosetta_STD_Disrupt, yerr=norm_STD_Disrupt, color='black', capsize=6, linestyle = '')
    #ax.errorbar(PyRosetta_WT_like, norm_WT_like, xerr=PyRosetta_STD_WT_like, yerr=norm_STD_WT_like, color='black', capsize=6, linestyle = '')

    '''
    for x, y, n in zip(PyRosetta_Disrupt, norm_Disrupt, N_Disrupt):
        if n=='V1414W' or n=='G1413W':
            plt.annotate(n, xy=(x,y), xytext=(-17, 17),
                     textcoords='offset points',
                     ha="right", va='bottom',
                     bbox=dict(boxstyle='round,pad=0.5',
                               fc='red', alpha=0.5),
                               arrowprops=dict(arrowstyle = '->',connectionstyle='arc3,rad=0',color='black'))
        else:
            plt.annotate(n, xy=(x,y), xytext=(57, -17),
                     textcoords='offset points',
                     ha="right", va='bottom',
                     bbox=dict(boxstyle='round,pad=0.5',
                               fc='red', alpha=0.5),
                               arrowprops=dict(arrowstyle = '->',connectionstyle='arc3,rad=0',color='black'))

    for x, y, n in zip(PyRosetta_WT_like, norm_WT_like, N_WT_like):
        if n=='T1416M' or n=='P1411M':
            plt.annotate(n, xy=(x,y), xytext=(57, 37),
                     textcoords='offset points',
                     ha="right", va='bottom',
                     bbox=dict(boxstyle='round,pad=0.5',
                               fc='blue', alpha=0.5),
                               arrowprops=dict(arrowstyle = '->',connectionstyle='arc3,rad=0',color='black'))
        elif n=='F1417Y' or n=='S1412M' or n=='G1413M':
            plt.annotate(n, xy=(x,y), xytext=(-17, -17),
                     textcoords='offset points',
                     ha="right", va='bottom',
                     bbox=dict(boxstyle='round,pad=0.5',
                               fc='blue', alpha=0.5),
                               arrowprops=dict(arrowstyle = '->',connectionstyle='arc3,rad=0',color='black'))
        else:
            plt.annotate(n, xy=(x,y), xytext=(57, -27),
                     textcoords='offset points',
                     ha="right", va='bottom',
                     bbox=dict(boxstyle='round,pad=0.5',
                               fc='blue', alpha=0.5),
                               arrowprops=dict(arrowstyle = '->',connectionstyle='arc3,rad=0',color='black'))
    '''
    '''
    for x, y, n in zip(PyRosetta_Disrupt, abs_pull_Disrupt, N_Disrupt):
        if n=='V1414W':
            plt.annotate(n, xy=(x,y), xytext=(-17, 17),
                     textcoords='offset points',
                     ha="right", va='bottom',
                     bbox=dict(boxstyle='round,pad=0.5',
                               fc='red', alpha=0.5),
                               arrowprops=dict(arrowstyle = '->',connectionstyle='arc3,rad=0',color='black'))
        else:
            plt.annotate(n, xy=(x,y), xytext=(57, -17),
                     textcoords='offset points',
                     ha="right", va='bottom',
                     bbox=dict(boxstyle='round,pad=0.5',
                               fc='red', alpha=0.5),
                               arrowprops=dict(arrowstyle = '->',connectionstyle='arc3,rad=0',color='black'))

    for x, y, n in zip(PyRosetta_WT_like, abs_pull_WT_like, N_WT_like):
        if n=='G1418A':
            plt.annotate(n, xy=(x,y), xytext=(57, 37),
                     textcoords='offset points',
                     ha="right", va='bottom',
                     bbox=dict(boxstyle='round,pad=0.5',
                               fc='blue', alpha=0.5),
                               arrowprops=dict(arrowstyle = '->',connectionstyle='arc3,rad=0',color='black'))
        elif n=='F1417Y' or n=='S1412M':
            plt.annotate(n, xy=(x,y), xytext=(-17, -17),
                     textcoords='offset points',
                     ha="right", va='bottom',
                     bbox=dict(boxstyle='round,pad=0.5',
                               fc='blue', alpha=0.5),
                               arrowprops=dict(arrowstyle = '->',connectionstyle='arc3,rad=0',color='black'))
        else:
            plt.annotate(n, xy=(x,y), xytext=(57, -27),
                     textcoords='offset points',
                     ha="right", va='bottom',
                     bbox=dict(boxstyle='round,pad=0.5',
                               fc='blue', alpha=0.5),
                               arrowprops=dict(arrowstyle = '->',connectionstyle='arc3,rad=0',color='black'))
    '''
    
    for x, y, n in zip(PyRosetta_Disrupt, FoldX_Disrupt, N_Disrupt):
            plt.annotate(n, xy=(x,y), xytext=(-17, 17),
                         textcoords='offset points',
                         ha="right", va='bottom',
                         bbox=dict(boxstyle='round,pad=0.5',
                                   fc='red', alpha=0.5),
                                   arrowprops=dict(arrowstyle = '->',connectionstyle='arc3,rad=0',color='black'))
    for x, y, n in zip(PyRosetta_WT_like, FoldX_WT_like, N_WT_like):
            plt.annotate(n, xy=(x,y), xytext=(-17, 17),
                         textcoords='offset points',
                         ha="right", va='bottom',
                         bbox=dict(boxstyle='round,pad=0.5',
                                   fc='blue', alpha=0.5),
                                   arrowprops=dict(arrowstyle = '->',connectionstyle='arc3,rad=0',color='black'))
    
    handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in ['red','blue']]
    labels= ['Disrupt', 'WT-like']
    plt.legend(handles, labels,prop={'size': 12})  
    
 
    
def R_Square(x_axis, y_axis):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_axis, y_axis)
    print ('slope ', slope)
    print ('intercept ', intercept)
    print ('r_value ', r_value)
    print ('p_value ', p_value)
    return (r_value,p_value)    
main()