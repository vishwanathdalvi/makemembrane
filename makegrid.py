# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 19:45:33 2015

@author: Ashwin
"""
import scipy
import scipy.linalg
import matplotlib.pyplot as plt

def getneighbours(pos):
    #Neighbours of nodes
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
    pos1 = (i-1, j+2)
    pos2 = (i+1, j-2)
    pos3 = (i  , j+1)
    pos4 = (i  , j-1)
    pos5 = (i-1, j+1)
    pos6 = (i+1, j-1)
    return set([pos1, pos2, pos3, pos4, pos5, pos6])

    
def getcellbonds(pos):
    [i, j] = pos
    bond1 = frozenset([(i  ,j  ,0),(i  ,j  ,1)])
    bond2 = frozenset([(i-1,j+2,0),(i-1,j+2,1)])
    bond3 = frozenset([(i  ,j  ,0),(i-1,j+1,1)])
    bond4 = frozenset([(i  ,j  ,1),(i  ,j+1,0)])
    bond5 = frozenset([(i-1,j+2,0),(i-1,j+1,1)])
    bond6 = frozenset([(i-1,j+2,1),(i  ,j+1,0)])
    return {bond1, bond2, bond3, bond4, bond5, bond6}
    
def getedge(cell1, cell2):
    #get the edge shared by the two cells
    [i1, j1] = cell1
    [i2, j2] = cell2
    if   (i1-i2)*(j1-j2) == -1:
        i1 = min([i1, i2])
        j1 = max([j1, j2])        
        node1 = (i1  , j1  , 1)
        node2 = (i1  , j1+1, 0)
    elif (i1 == i2) and (abs(j1-j2) == 1):
        i1 = i1-1
        j1 = min([j1, j2])+2
        node1 = (i1  , j1  , 1)
        node2 = (i1+1, j1-1, 0)
    elif (i1-i2)*(j1-j2) == -2:
        i1 = min([i1, i2])
        j1 = max([j1, j2])
        node1 = (i1  , j1  , 0)
        node2 = (i1  , j1  , 1)
    return {node1, node2}
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

def getcells(setofnodes):
    setofporecells = set([])
    maxx = 0.0
    maxy = 0.0
    for node in setofnodes:
        [i,j,k] = node
        cell = (i,j)
        setofporecells.add(cell)
        x, y = getcellposition(cell)
        maxx = max([abs(x), maxx])
        maxy = max([abs(y), maxy])
    return setofporecells, maxx, maxy
    
def isclose(position, cell):
    x, y = getcellposition(cell)
    px, py = position
    dist = scipy.linalg.norm([x-px, y-py])
    if dist <= scipy.sqrt(3.0)*0.5:
        boolfound = True
    else:
        boolfound = False
    return boolfound

def getindex(position):
    [px, py] = position
    cell = (0,0)
    setoftraversedcells = set([])
    setofbordercells = set([cell])
    boolfound = False
    while not boolfound:
        newsetofbordercells = setofbordercells.union(set([]))
        for cell in setofbordercells:
            boolfound = isclose(position, cell)
            if boolfound:
                break
            else:            
                setoftraversedcells.add(cell)
                newsetofbordercells.remove(cell)
                neighbours = getcellneighbours(cell)
                newsetofbordercells = newsetofbordercells.union(neighbours)
        setofbordercells = newsetofbordercells.difference(setoftraversedcells)
    return cell

def getporeborderbonds(pore):
    setofborderbonds, setofinnerbonds, setofbordernodes, setofporenodes = set([]), set([]), set([]), set([])
    setoftraversedbonds = set([])
    for cell in pore:
        setofcellbonds = getcellbonds(cell).difference(setoftraversedbonds)
        for bond in setofcellbonds:
            setoftraversedbonds.add(bond)
            setofedgeneighbours = getedgeneighbours(bond).difference(pore)
            if len(setofedgeneighbours) == 1:
                setofborderbonds.add(bond)
                [node1, node2] = bond
                setofbordernodes.add(node1)
                setofbordernodes.add(node2)
            else:
                [node1, node2] = bond                
                setofporenodes.add(node1)
                setofporenodes.add(node2)
            if len(setofedgeneighbours) == 0:
                setofinnerbonds.add(bond)
    setofporenodes = setofporenodes.difference(setofbordernodes)
    return setofborderbonds, setofinnerbonds, setofbordernodes, setofporenodes
            
#def getborderpositions(setofborderbonds):
         
    
    

def plotpores(ax, listofpores, setofbonds):
    for pore in listofpores:
        porepos = [getcellposition(p) for p in pore]
        x = [p[0] for p in porepos]
        y = [p[1] for p in porepos]
        ax.plot(x, y, 'ro')
        for cell in pore:
            setofcellbonds = getcellbonds(cell).intersection(setofbonds)
            for bond in setofcellbonds:
                node1, node2 = bond
                [x1, y1] = getposition(node1)
                [x2, y2] = getposition(node2)
                ax.plot([x1, x2],[y1, y2], 'k')
        setofborderbonds, setofinnerbonds, setofbordernodes, setofporenodes = getporeborderbonds(pore)
        for bond in setofborderbonds:
            node1, node2 = bond
            [x1, y1] = getposition(node1)
            [x2, y2] = getposition(node2)
            ax.plot([x1, x2],[y1, y2], color = 'blue', linewidth = 2.0, alpha = 0.5)            
        for bond in setofinnerbonds:
            node1, node2 = bond
            [x1, y1] = getposition(node1)
            [x2, y2] = getposition(node2)
            ax.plot([x1, x2],[y1, y2], color = 'cyan', linewidth = 2.0, alpha = 0.5)
        for node in setofbordernodes:
            [x, y] = getposition(node)
            ax.plot([x],[y],'bs')
        for node in setofporenodes:
            [x, y] = getposition(node)
            ax.plot([x],[y], 'bo')
        
def getopenneighbours(cell, setofbonds, maxx, maxy):
    neighbours = getcellneighbours(cell)
    openneighbours = set([])
    boolborder = False
    for neighbour in neighbours:
        x, y = getcellposition(neighbour)
        if abs(x) > maxx or abs(y) > maxy:
            boolborder = True
        bond = getedge(cell, neighbour)
        if bond not in setofbonds:
            openneighbours.add(neighbour)
    return openneighbours, boolborder

    
def getpores(setofnodes, setofbonds):
    setofcells, maxx, maxy = getcells(setofnodes)
    listofpores = [] 
    setoftraversedcells = set([])
    for cell in setofcells:
        if cell not in setoftraversedcells:
            boolstop = False
            pore = set([cell])
            setofbordercells = set([cell])
            while not boolstop:
                newsetofbordercells = setofbordercells.union(set([]))
                for bordercell in setofbordercells:
                    openneighbours, boolborder = getopenneighbours(bordercell, setofbonds,maxx, maxy)
                    if not boolstop:
                        pore.add(bordercell)
                        newsetofbordercells = newsetofbordercells.union(openneighbours).difference(pore)
                        
                    if boolborder:
                        boolstop = True
                        break
                setofbordercells = newsetofbordercells.difference(pore)
                if len(setofbordercells) == 0:
                    boolstop = True
            setoftraversedcells = setoftraversedcells.union(pore)     
            if not boolborder:
                listofpores.append(pore)          
    return listofpores

