# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 19:45:33 2015

@author: Ashwin
"""
import scipy
import matplotlib.pyplot as plt

def getneighbours(pos):
    [i, j, k] = pos  #Left = 0 and Right = 1
    kn = 0 if k == 1 else 1
    if kn == 0:
        pos1 = (i  , j  , kn)
        pos2 = (i  , j+1, kn)
        pos3 = (i+1, j-1, kn)
    elif kn == 1:
        pos1 = (i  , j  , kn)
        pos2 = (i  , j-1, kn)
        pos3 = (i-1, j+1, kn)
    return [pos1, pos2, pos3]

def getposition(pos):
    [i, j, k] = pos
    px = 3*(i + 0.5*j) - 0.5*(-2*k+1)
    py = scipy.sqrt(3)*0.5*j
    return [px, py]
    





        
    