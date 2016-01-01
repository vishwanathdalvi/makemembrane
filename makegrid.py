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

def getcellneighbours(pos):
    [i, j] = pos 
    pos1 = [i-1, j+2]
    pos2 = [i+1, j-2]
    pos3 = [i  , j+1]
    pos4 = [i  , j-1]
    pos5 = [i-1, j+1]
    pos6 = [i+1, j-1]
    return set([pos1, pos2, pos3, pos4, pos5, pos6])

def getcellbonds(pos):
    [i, j] = pos
    bond1 = {(i  ,j  ,0),(i  ,j  ,1)}
    bond2 = {(i-1,j+2,0),(i-1,j+2,1)}
    bond3 = {(i  ,j  ,0),(i-1,j+1,1)}
    bond4 = {(i  ,j  ,1),(i  ,j+1,0)}
    bond5 = {(i-1,j+2,0),(i-1,j+1,1)}
    bond6 = {(i-1,j+2,1),(i  ,j+1,0)}
    return [bond1, bond2, bond3, bond4, bond5, bond6]
    
def getcells(setofnodes):
    setofporecells = set([])
    for node in setofnodes:
        [i,j,k] = node
        cell = (i,j)
        setofporecells.add(cell)
    return setofporecells

def getedgeneighbours(bond):
    #Get the two cells sharing this bond
    [node1, node2] = bond
    [i_1,j_1,k1] = node1
    [i_2,j_2,k2] = node2
    if k1 == 1:
        i1, j1 = i_1, j_1
        i2, j2 = i_2, j_2
    else:
        i1, j1 = i_2, j_2
        i2, j2 = i_1, j_1
    if   (i1,j1) == (i2  ,j2  ):
        cell1, cell2 = (i1  , j1  ), (i1+1, j1-2)
    elif (i1,j1) == (i2  ,j2-1):
        cell1, cell2 = (i1  , j1  ), (i1+1, j1-1)
    elif (i1,j1) == (i2-1,j2+1):
        cell1, cell2 = (i1+1, j1-1), (i1+1, j1-2)
    return {cell1, cell2}
        
def getcellposition(cell):
    [i, j] = cell
    px = 3*(i+0.5*j)
    py = 0.5*3**0.5*(j+1)
    return px, py



        
    