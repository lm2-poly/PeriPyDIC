# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import csv
import numpy as np

## Class handels the output for CSV files
class OutputCSV():
     ## Constructor 
     # @param outType Type of the output 
     # @param dataType The data, which is wriiten to the CSV file 
     # @param inputFile The file where the output is written
     def __init__(self,outType,dataType,inputFile):
         ## Type of the output
         self.outType = outType
         ## Filename for the output file
         self.inputFile = inputFile
         ## Type of the written data
         self.dataType = dataType
         
     ## Writes the data to the CSV file 
     # @param PD_deck The object containing the configuration of the yaml file
     # @param problem The object containing the computational vlaues
     def write(self, PD_deck,problem):
        with open(self.inputFile, 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=' ')
            header = []
            
            if self.dataType == "Position":
                
                header.append("#Time")
                for x_i in range(0,PD_deck.num_nodes_x):
                    header.append("Id"+str(x_i))        
                spamwriter.writerow(header)
            
                for t_n in range(0, PD_deck.time_steps):
                    spamwriter.writerow( np.insert(problem.y[:,t_n] , 0,  t_n*PD_deck.delta_t   ))
        
