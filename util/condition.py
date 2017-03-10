import os
import sys
import csv
import numpy as np


class ConditionFromFile():            

    boundary_volume = 0.0
    
    def __init__(self,type,inFile,value,volume):
        self.id = self.readCondition(inFile,volume)
        self.type = type
        self.value = float(value)
        self.force_density = self.value / self.boundary_volume
                        
    def readCondition(self, inFile,volume):
        """ Reads the ids from the inFile where this condition should be 
            applied.
            
            Args:
                inFile: File name of the CSV file with the ids of the nodes
                    for this boundary condition
            volume: Volume of the nodes
            
            Returns:
                id: The ids read from the inFile
        """
        if not os.path.exists(inFile):
            print "Error: Could not find " + inFile
            sys.exit(1)

        with open(inFile, 'r') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=' ')
            # Skip the first line, because is the header
            next(spamreader)
            length = len(list(spamreader))
            id = np.empty(length)
            csvfile.seek(0)
            next(spamreader)
            i = 0
            for row in spamreader:
                id[i] = int(row[0])
                self.boundary_volume += volume[int(id[i])]
                i += 1
                
            return id
      
