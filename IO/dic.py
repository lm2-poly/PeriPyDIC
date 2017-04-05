# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca

import csv
import numpy as np

class DICreader2D():
    
    """
    A class for reading an providing files from 2 dimensional dic
    """

    def __init__(self, path):
        self.dim = 2
        self.data = []
        
        self.read(path)
        self.sortData()
        self.extractData()
        self.determineUnitHorizon()
        

    def read(self, path):
        """
        Read the values provided by the dic and stores it to the data array
        `path`  Path and appended file name for the csv file to proceed
        """
        self.length = 0
        with open(path, 'rb') as csvfile:
            csvreader = csv.reader(csvfile, delimiter=',')
            next(csvreader, None)
            for row in csvreader:
                self.data.append(np.array(map(float, row)))
                self.length += 1
                
    def sortData(self):
        """
        Sorts the data from the csv file with respect to the first and second pixel
        """
        self.data.sort(key=lambda x: (x[0], x[1]))
        
        
    def determineUnitHorizon(self):
        # Find unique values for x:
        sorted_x_set = sorted(set(self.x))
        sorted_y_set = sorted(set(self.y))
        
        self.delta_x = abs(sorted_x_set[1] - sorted_x_set[0])
        
    def extractData(self):
        self.x = np.zeros((self.length))
        self.dx = np.zeros((self.length))
        self.y = np.zeros((self.length))
        self.dy = np.zeros((self.length))
        
        for i in range(0, len(self.data)):
            self.x[i] = self.data[i][0]
            self.y[i] = self.data[i][1]
            
            self.dx[i] = self.data[i][3]
            self.dy[i] = self.data[i][4]
            
        del self.data
        
        self.nodes = np.empty((self.length, self.dim))
        self.nodes[:,0] = self.x
        self.nodes[:,1] = self.y
        
