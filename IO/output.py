# -*- coding: utf-8 -*-
"""
Created on Thu Mar 09 08:51:08 2017

@author: Rolland
"""
import csv
import numpy as np

class OutputCSV():
     
     def __init__(self,outType,dataType,inputFile):
         self.outType = outType
         self.inputFile = inputFile
         self.dataType = dataType
         

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
        