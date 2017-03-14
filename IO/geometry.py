#-*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import numpy as np
import csv
import os
import sys

## Class handeling the discrete nodes
class Geometry():
    ## Constructor
    def __init__(self):
        ## Positions in x-direction
        self.pos_x = np.array(0)
        ## Positions in y-direction
        self.pos_y = np.array(0)
        ## Positions in z-direction
        self.pos_z = np.array(0)
        ## Volume per node
        self.volumes = np.array(0)
        ## Density per node
        self.density = np.array(0)
        ## Amount of nodes
        self.amount = 0
        
        self.volume_boundary = 0.
    
    ## Read the positions, volume, and density of the nodes from the inFile.
    # @param dim Dimension of the nodes
    # @param inFile CSV file with the geometry
    def readNodes(self,dim,inFile):
        if not os.path.exists(inFile):
                print "Error: Could not find " + inFile
                sys.exit(1)
        #Dimension of the problem
        self.dim = dim
        with open(inFile, 'r') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=' ')
            #Skip the first line, because is the header
            next(spamreader)
            length = len(list(spamreader))
            csvfile.seek(0)
            next(spamreader)
            
            if dim <= 1:
                self.pos_x = np.empty(length)
            if dim <= 2:
                self.pos_y = np.empty(length)
            if dim <= 3:
                self.pos_z = np.empty(length)
            
            self.volumes = np.empty(length)
            i = 0
            for row in spamreader:
                if dim >= 1:
                    self.pos_x[i] = float(row[1])
                if dim >= 2:
                    self.pos_y[i] = float(row[2])
                if dim == 3:
                    self.pos_z[i] = float(row[3])
                
                self.volumes[i] = float(row[dim+1])
                i +=1
        self.generateNodes()
    ## Computes the min distance between all nodes
    # @param direction
    # @return Minimal direction    
    def getMinDist(self,direction):
        tmp = float('inf')
        if direction == 1:
            for i in range(0,len(self.nodes)):
                for j in range(0,len(self.nodes)):
                    if i != j:
                        val = abs(self.nodes[j]-self.nodes[i])
                        if val < tmp:
                            tmp = val
        if direction == 2:
            for i in range(0,len(self.nodes[1])):
                for j in range(0,len(self.nodes[1])):
                    if i != j:
                        val = abs(self.nodes[1][j]-self.nodes[1][i])
                        if val < tmp:
                            tmp = val
        if direction == 3:
            for i in range(0,len(self.nodes[1])):
                for j in range(0,len(self.nodes[1])):
                    if i != j:
                        val = abs(self.nodes[2][j]-self.nodes[2][i])
                        if val < tmp:
                            tmp = val
        return tmp
                    
    def generateNodes(self):
           
        if self.dim == 1:
            self.nodes = np.array((self.pos_x),dtype=np.double)
            del self.pos_x
            self.amount = len(self.nodes)
        if self.dim == 2:
            self.nodes = np.array((self.pos_x,self.pos_y),dtype=np.double)
            del self.pos_x
            del self.pos_y
            self.amount = len(self.nodes[0])
        if self.dim == 3:
            self.nodes = np.array((self.pos_x,self.pos_y,self.pos_z),dtype=np.double)
            del self.pos_x
            del self.pos_y
            del self.pos_z
            self.amount = len(self.nodes[0])
            
