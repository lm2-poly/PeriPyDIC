# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 12:50:28 2015

@author: ilyass.tabiai@gmail.com
@author: rolland.delorme@gmail.com
@author: diehl@ins.uni-bonn.de
"""
import csv
import numpy as np

class dic_reader_2D():
  """
  A class for reading an providing files from 2 dimensional dic
  """
 
  def __init__(self,path):
    self.data = []
    self.read(path)
    self.sortdata()
  
  def read(self,path):
      """
      Read the values provided by the dic and sotres it to the data array
      
      `path`  Path and appended file name for the csv file to proced
      """
      with open(path,'rb') as csvfile:
	csvreader = csv.reader(csvfile, delimiter=',')
	next(csvreader,None)
	for row in csvreader:
	   self.data.append(np.array(map(float, row)))
	    
  def sortdata(self):
    """
    Sorts the data from the csv file with respect to the first and second pixel
    """
    self.data.sort(key = lambda x: (x[13], x[14]))
    print len(self.data)
    for i in range(0,len(self.data)):
      print self.data[i]
    
  
  














