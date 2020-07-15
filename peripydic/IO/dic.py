# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrickdiehl@lsu.edu
import csv
import numpy as np
from . import deck
import matplotlib.pyplot as plt

## A class for reading VIC3D CSV grid exports and converting them into list
# objects usable by other classes
class DICreader2D():

    ## Constructor
    # Reads the CSV file provided, determines the unit horizon of this data and
    # and creates node initial and node actual attributes of the class usable
    # by the other modules
    # @param deck Deck object containing input data from the .yaml file
    def __init__(self, deck):
        ## Dimension of the problem (2D for DIC)
        self.dim = deck.dim
        ## Temporary variable internal to this class
        self.data = []
        self.sigma_column = deck.sigma_column

        self.read(deck.filepath + "/" + deck.filename)
        self.extractData(deck)
        self.determineUnitHorizon()

    ## Read the values provided by the dic and stores it to the data array
    #`path`  Path and appended file name for the csv file to proceed
    # @param path The path to the VIC3D CSV export data
    def read(self, path):
        ## Total number of points in the DIC CSV file
        self.length = 0
        with open(path, 'r') as csvfile:
            csvreader = csv.reader(csvfile, delimiter=',')
            next(csvreader, None)
            for row in csvreader:
                self.data.append(np.array(list(map(float, row))))
                self.length += 1

    ## Find unique values for x
    # @param self Object pointer
    def determineUnitHorizon(self):
        self.delta_x = np.sqrt(np.power(self.x[1]-self.x[0],2)+np.power(self.y[1]-self.y[0],2))

    ## Stores the data extracted from the CSV file in objects which can be
    # manipulated by other modules
    # @param self Object pointer
    def extractData(self,deck):
        ## Vector containing x nodes from DIC
        self.x = np.zeros((self.length))
        dx = np.zeros((self.length))
        ## Vector containing y nodes from DIC
        self.y = np.zeros((self.length))
        dy = np.zeros((self.length))
        ## Vector containing columes for each point
        self.volumes = np.zeros((self.length))
        ## Strain from DIC
        self.strain = np.zeros((self.length,3))

        if deck.filetype == "VIC3D":
            for i in range(0, len(self.data)):
                # Remove values for which the confidence is -1.0
                if self.data[i][self.sigma_column] == -1.:
                    pass
                else:

                    self.x[i] = self.data[i][0]
                    if self.dim == 2:
                        self.y[i] = self.data[i][1]

                    dx[i] = self.data[i][3]
                    if self.dim == 2:
                        dy[i] = self.data[i][4]

                    self.strain[i][0] = self.data[i][6]
                    self.strain[i][1] = self.data[i][7]
                    self.strain[i][2] = self.data[i][8]

                    self.volumes[i] = deck.dic_volume * deck.thickness

            del self.data
            ## Nodes initial positions
            self.nodes = np.empty((self.length, self.dim))
            self.nodes[:,0] = self.x 
            if self.dim == 2:
                self.nodes[:,1] = self.y 
            ## Nodes actual positions
            self.act = np.empty((self.length, self.dim))
            self.act[:,0] = self.x  + dx
            if self.dim == 2:
                self.act[:,1] = self.y + dy


        if deck.filetype == "mudic":

            for i in range(0, len(self.data)):
                self.x[i] = self.data[i][1]
                if self.dim == 2:
                    self.y[i] = self.data[i][2]
                
                dx[i] = self.data[i][4]

                if self.dim == 2:
                    dy[i] = self.data[i][5]
                    
                self.volumes[i] = self.data[i][3] * deck.thickness

                self.strain[i][0] = self.data[i][8]
                self.strain[i][1] = self.data[i][9]
                self.strain[i][2] = self.data[i][10]


            del self.data
            ## Nodes initial positions
            self.nodes = np.empty((self.length, self.dim))
            self.nodes[:,0] = self.x -dx
            if self.dim == 2:
                self.nodes[:,1] = self.y -dy
            ## Nodes actual positions
            self.act = np.empty((self.length, self.dim))
            self.act[:,0] = self.x  #+ dx
            if self.dim == 2:
                self.act[:,1] = self.y #+ dy

