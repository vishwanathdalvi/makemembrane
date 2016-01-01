# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 20:11:14 2015

@author: Ashwin
"""

import makegrid as mkg
import scipy
import matplotlib.pyplot as plt
import random
from random import random as rnd
from random import randint as rint
from random import shuffle, sample

pause = 0.001
'''
C   Free Chloride
CB  Singly bound chloride
C2B Doubly bound chloride
C3B Triply (fully) bound chloride
A   Free Amine
AB  Singly bound amine
A2B Doubly (fully) bound amine
P   Protonated singly bound amine 
'''
dict_symbols = {'C' :'yo', 
                'A' :'y^',
                'CB' :'mo',
                'C2B':'mo',
                'C3B':'go',
                'AB' :'r^',
                'A2B':'c^',
                'P'  :'r^'}

reactionsets = [{'C'  ,'A' },
                {'C'  ,'AB'},
                {'CB' ,'A' },
                {'CB' ,'AB'},
                {'C2B','A' },
                {'C2B','AB'},
                {'P'  ,'A' }]

def isreactionold(typ1, typ2):
    if typ1 in ['C3B','A2B']:
        return False
        
    if typ2 in ['C3B','A2B']:
        return False
        
    if typ1 in ['C','CB','C2B']:
        if typ2 in ['C','CB','C2B','P']:
            return False
        if typ2 in ['A','AB']:
            return True
    if typ2 in ['C','CB','C2B']:
        if typ1 in ['C','CB','C2B','P']:
            return False
        if typ1 in ['A','AB']:
            return True

    if typ1 in ['A','AB']:
        if typ2 in ['C', 'CB', 'C2B', 'P']:
            return True
        if typ2 in ['A','AB']:
            return False
    if typ2 in ['A','AB']:
        if typ1 in ['C', 'CB', 'C2B', 'P']:
            return True
        if typ1 in ['A','AB']:
            return False
            
    
    if typ1 == 'P':
        if typ2 == 'A':
            return True
        else:
            return False
    if typ2 == 'P':
        if typ1 == 'A':
            return True
        else:
            return False

def isreaction(typ1, typ2):
    sett = {typ1, typ2}
    if sett in reactionsets:
        return True
    else:
        return False
    
    

dict_transform = {'C'  :'CB'  ,
                  'CB' :'C2B' ,
                  'C2B':'C3B' ,
                  'A'  :'P'   ,
                  'AB' :'A2B' ,
                  'P'  :'AB'  }

typUNreactivenodes = ['A2B','C3B']


def plot(ax, setofnodes, listofbonds, dictattributes, bool_test = False):
    ax.cla()
    for bond in listofbonds:
        [node1, node2] = bond
        [x1, y1] = mkg.getposition(node1)
        [x2, y2] = mkg.getposition(node2)
        ax.plot([x1, x2], [y1, y2], color = 'black')
    for node in setofnodes:
        attr = dictattributes[node]
        [x, y] = mkg.getposition(node)
        ax.plot([x], [y], dict_symbols[attr[0]])
        
        if bool_test:
            for neighbour in attr[1]:
                [xn, yn] = mkg.getposition(neighbour)
                ax.plot([x, xn], [y, yn], color = 'black')

def plotpores(ax, listofpores):
    for pore in listofpores:
        porepos = [mkg.getcellposition(p) for p in pore]
        x = [p[0] for p in porepos]
        y = [p[1] for p in porepos]
        ax.plot(x, y, 'ro')
        ax.plot(x, y, 'r--')

def getopenborders(cell, listofbonds):
    cellborders = mkg.getcellbonds(cell) #get the bonds surrounding that cell
    openborders = [bond for bond in cellborders if bond not in listofbonds]#Select only those borders that don't have bonds 
    return openborders
    
def tracepore(pore, nopbor, setofnodes, listofbonds, setofcells, setoftraversedcells):
    print nopbor
    cell = pore[-1]
    openborders = getopenborders(cell, listofbonds)
    for bond in openborders:
        edgecellpair = mkg.getedgeneighbours(bond)
        cellnew = edgecellpair.difference(set(pore))
        if len(cellnew) == 1:
            cellnew = list(cellnew)[0]
            if cellnew not in setoftraversedcells:
                pore.append(cellnew)
                setoftraversedcells.add(cellnew)
                openborders = getopenborders(cellnew, listofbonds)
                nopbor = max([nopbor, len(openborders)])
                if len(openborders) < 6:
                    pore, nopbor = tracepore(pore, nopbor, setofnodes, listofbonds, setofcells, setoftraversedcells)
    return pore, nopbor
    
def getpores(setofnodes, listofbonds):
    setofcells = mkg.getcells(setofnodes)
    listofpores = [] 
    setoftraversedcells = set([])
    for cell in setofcells:
        setoftraversedcells.add(cell)
        openborders = getopenborders(cell, listofbonds)
        nopbor = len(openborders)
        pore, nopbor = tracepore([cell], nopbor, setofnodes, listofbonds, setofcells, setoftraversedcells)
        cell = pore[-1]
        openborders = getopenborders(cell, listofbonds)
        if nopbor < 6:   
            listofpores.append(pore)
    return listofpores
        
            
class Simulation:
    def __init__(self):
        pos0 = (0,0,0)
        self.setofnodes = {pos0}
        self.listofbonds = []
        
        self.setofreactivenodes = {pos0}
        self.dict_attributes = {pos0:['C',set([])]} #First attribute is type, second is set of bonded neighbours
        self.setofbordernodes = set(mkg.getneighbours(pos0))
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
    def addreactant(self, p, boolplot=False): #p is probability of getting an amine
        r = rnd()
        if r < p:
            moltyp =  'A'
        else:
            moltyp = 'C'
        node = sample(self.setofbordernodes, 1)[0]
        if boolplot:
            self.ax.clear()
            plot(self.ax, self.setofnodes, self.dict_attributes)
            [px, py] = mkg.getposition(node)
            if moltyp == 'A':
                self.ax.plot([px],[py],'y^')
            elif moltyp == 'C':
                self.ax.plot([px],[py],'yo')
            self.fig.canvas.draw()
            plt.pause(pause)
        return moltyp, node 
    
    def internalreact(self):
        boolreaction = True
        while boolreaction:
            boolreaction = False
            for reactnode in self.setofreactivenodes:
                neighbours = set(mkg.getneighbours(reactnode))
                reactiveneighbours = self.setofreactivenodes.intersection(neighbours)
                typ = self.dict_attributes[reactnode][0]
                for reactiveneighbour in reactiveneighbours:
                    if reactiveneighbour not in self.dict_attributes[reactnode][1]: #Make sure neighbour is not already bonded
                        typneighbour = self.dict_attributes[reactiveneighbour][0]
                        boolreact = isreaction(typ,typneighbour)
                        if boolreact:
                            boolreaction = True
                            
                            typtrans = dict_transform[typ]
                            self.dict_attributes[reactnode][0] = typtrans
                            self.dict_attributes[reactnode][1].add(reactiveneighbour)
    
                            typneighbourtrans = dict_transform[typneighbour]
                            self.dict_attributes[reactiveneighbour][0] = typneighbourtrans
                            self.dict_attributes[reactiveneighbour][1].add(reactnode)
    
                            if typtrans in typUNreactivenodes:
                                self.setofreactivenodes.remove(reactnode)
    
                            if typneighbourtrans in typUNreactivenodes:
                                self.setofreactivenodes.remove(reactiveneighbour)
                            
                            bond = {reactnode, reactiveneighbour}
                            if bond not in self.listofbonds:
                                self.listofbonds.append(bond)
    
                            break
                if boolreaction:
                    break
    def react(self, p, boolplot=False):                
        
        moltyp, node = self.addreactant(p, boolplot = boolplot)
        setofneighbours = set(mkg.getneighbours(node))
        setofreactiveneighbours = list(self.setofreactivenodes.intersection(setofneighbours))
        shuffle(setofreactiveneighbours)
        
        booladd = False
        for reactnode in setofreactiveneighbours:
            typ = self.dict_attributes[reactnode][0]
            boolreact = isreaction(moltyp,typ)

            if boolreact:
                booladd = True
                newtyp = dict_transform[typ]
                self.dict_attributes[reactnode][0] = newtyp
                if newtyp in typUNreactivenodes:
                    self.setofreactivenodes.remove(reactnode)
                if typ == 'P':
                    booladd = False
                break

        if booladd:
            self.dict_attributes[reactnode][1].add(node) #New node added to list of bonded neighbours of reacted node
            bond = {reactnode, node}
            if bond not in self.listofbonds:
                self.listofbonds.append(bond)
                
            self.setofbordernodes.remove(node) #Border node removed from set of border nodes
            self.setofnodes.add(node) #Node added to set of nodes
            transmoltyp = dict_transform[moltyp]
            self.dict_attributes[node] = [transmoltyp, {reactnode}] #reactnode added to list of bonded neighbours
            
            if transmoltyp not in typUNreactivenodes:
                self.setofreactivenodes.add(node)
            
            bordernodes = set(mkg.getneighbours(node))
            invalidbordernodes = self.setofnodes.intersection(bordernodes)
            newbordernodes = bordernodes.difference(invalidbordernodes)
            self.setofbordernodes = self.setofbordernodes.union(newbordernodes)

        if boolplot:
            self.ax.clear()
            plot(self.ax, self.setofnodes, self.dict_attributes)
            self.fig.canvas.draw()
            plt.pause(pause)
    
    def cure(self): #Convert protonated amines into reactive ones
        for node in self.setofnodes:
            attr = self.dict_attributes[node]
            typ = attr[0]
            if typ == 'P':
                attr[0] = 'AB'
        self.internalreact()
    
    def simulate(self, p, n=10, boolplot = False):
        for i in xrange(n):
            self.internalreact()
            self.react(p)
        self.cure()
        self.listofpores = getpores(self.setofnodes, self.listofbonds)
        if boolplot:
            plot(self.ax, self.setofnodes, self.listofbonds, self.dict_attributes)
            plotpores(self.ax, self.listofpores)
            self.fig.canvas.draw()
        
if __name__ == "__main__":
    sim = Simulation()


 
