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
 
  def __init__(self,path):
    self.data = []
    self.read(path)
    self.sortdata()
  
  def read(self,path):
      with open(path,'rb') as csvfile:
	csvreader = csv.reader(csvfile, delimiter=',')
	next(csvreader,None)
	for row in csvreader:
	   self.data.append(np.array(map(float, row)))
	    
  def sortX(self,item):
      return item[0]
      
  def sortY(self,item):
      return item[1]
  
  def sortZ(self,item):
      return item[2]
      
  def sortdata(self):
    self.data.sort(key = lambda x: (x[13], x[14]))
    print len(self.data)
    for i in range(0,len(self.data)):
      print self.data[i]
    
  
  














