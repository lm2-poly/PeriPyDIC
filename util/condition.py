#-*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import os
import sys
import csv
import numpy as np

## Class for storing the conditions from the yaml file
class ConditionFromFile():

    ## Volume of the boundary elements
    boundary_volume = 0.0

    ## Constructor
    # @param cType Type of the condition, e.g. Force or Displacement
    # @param inFile Path to the CSV file with the respective node ids
    # @param value The value, which is applied
    # @param volume The volume of the nodes
    # @param direction The direction where the conditions is applied
    def __init__(self,cType,inFile,value,volume,direction):
        ## Ids of the node where this condition is applied
        self.id = self.readCondition(inFile,volume)
        ## Type of the condition (Force or Displacement)
        self.type = cType
        ## Value in Newton or Meter
        self.value = float(value)
        ## Force density
        self.force_density = self.value / self.boundary_volume
        ## Direction of the condition
        self.direction = direction

    ##Reads the ids from the inFile where this condition should be applied.
    # @param inFile File name of the CSV file with the ids of the nodes
    # @param volume The volume of the nodes
    # @return The ids read from the inFile
    def readCondition(self, inFile,volume):
        if not os.path.exists(inFile):
            print "Error: Could not find " + inFile
            sys.exit(1)

        with open(inFile, 'r') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=' ')
            # Skip the first line, because is the header
            next(spamreader)
            length = len(list(spamreader))
            ids = np.empty(length)
            csvfile.seek(0)
            next(spamreader)
            i = 0
            for row in spamreader:
                ids[i] = int(row[0])
                self.boundary_volume += volume[int(ids[i])]
                i += 1

            return ids

