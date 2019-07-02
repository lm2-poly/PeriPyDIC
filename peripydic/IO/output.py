# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrickdiehl@lsu.edu
import csv
import numpy as np

## Class handles the output for CSV files
class OutputCSV():
    ## Constructor
    # @param outType Type of the output
    # @param dataType The data, which is written to the CSV file
    # @param inputFile The file where the output is written
    def __init__(self, outType, dataType, inputFile):
        ## Type of the output
        self.outType = outType
        ## Filename for the output file
        self.inputFile = inputFile
        ## Type of the written data
        self.dataType = dataType
    ## Writes the data to the CSV file
    # @param deck The object containing the configuration of the yaml file
    # @param problem The object containing the computational values
    def write(self, deck, problem):
        with open(self.inputFile, 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=' ', quotechar='|')
            header = []
            header.append("#Time")
            header.append("ID")
            if self.dataType == "Position":
                pos = ['X','Y','Z']
                for i in range(0, deck.dim):
                    header.append(pos[i])
                spamwriter.writerow(header)
                for t in range(0, deck.time_steps):
                    s = [t]
                    for i in range(0,deck.num_nodes):
                        s.append(i)
                        if deck.dim >= 1:
                            s.append(problem.y[i][0][t])
                        if deck.dim >= 2:
                            s.append(problem.y[i][1][t])
                        if deck.dim >= 3:
                            s.append(problem.y[i][2][t])
                        spamwriter.writerow(s)
                        s = [t]

