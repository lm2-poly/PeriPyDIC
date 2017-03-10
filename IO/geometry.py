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
       
        self.volume_boundary = 0.
    
    ## Read the positions, volume, and density of the nodes from the inFile.
    # @param dim Dimension of the nodes
    # @param inFile CSV file with the geometry
    def readNodes(self,dim,inFile):
        if not os.path.exists(inFile):
                print "Error: Could not find " + inFile
                sys.exit(1)
        
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
    
    ## Computes the min distance between all nodes
    # @param direction
    # @return Minimal direction    
    def getMinDist(self,direction):
        tmp = float('inf')
        if direction == 1:
            for i in range(0,len(self.pos_x)):
                for j in range(0,len(self.pos_x)):
                    if i != j:
                        val = abs(self.pos_x[j]-self.pos_x[i])
                        if val < tmp:
                            tmp = val
        if direction == 2:
            for i in range(0,len(self.pos_y)):
                for j in range(0,len(self.pos_y)):
                    if i != j:
                        val = abs(self.pos_y[j]-self.pos_y[i])
                        if val < tmp:
                            tmp = val
        if direction == 3:
            for i in range(0,len(self.pos_z)):
                for j in range(0,len(self.pos_z)):
                    if i != j:
                        val = abs(self.pos_z[j]-self.pos_z[i])
                        if val < tmp:
                            tmp = val
        return tmp
                    
        
