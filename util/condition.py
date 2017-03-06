import os
import sys
import csv
import numpy as np

class Condition():

    def __init__(self, data,volume_boundary):
        self.type = data["Type"]
        self.value = data["Value"]
        if "File" in data:
            self.id = self.readCondition(data["File"])
        self.symmetry = data["Symmetry"]
        if self.type == "Force":
            if self.symmetry == True:
                self.force_density = self.value / volume_boundary * 2.0
            else:
                self.force_density = self.value / volume_boundary

    def readCondition(self, inFile):
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
			i += 1
        return id
