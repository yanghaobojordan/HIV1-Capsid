# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 11:58:43 2020

@author: hyang
"""

def main():
    f=open("Folding_WT_8_dock_output.fasc","r")
    low=0
    structure=''
    f.readline()
    for i in f:
        i=i.rstrip().split()
        if float(i[3])<low:
            low=float(i[3])
            structure=i[1]
    print (low, structure)
main()