 #-*- coding: utf-8 -*-
"""
This class reads the geometries files and handles the specified
boundary conditions. In addition it could be use to discretize 
regular structures.
@author: patrick.diehl@polymtl.ca
""" 

import numpy as np

class Geometry():
    
    """
    Positions in x-direction
    """
    pos_x = np.array(0)
    
    """
    Positions in y-direction
    """
    pos_y = np.array(0)
    
    """
    Positions in z-direction
    """
    pos_z = np.array(0)
    
    """
    Volume per node
    """
    vol = np.array(0)
    
    """
    Density per node
    """
    density = np.array(0)
    
    volume_boundary = 0.
    
    def genrateGrid(dim,data,horizon_factor):
        if dim == 1:
                    generateGrid1D(data)
                    generateVolume1D(data,horizon_factor)
            
    
    def generateGrid1D():
        x = np.zeros(data["Nodes_X"])
        for i in range(0, data["Nodes_X"]):
            x[i] = (i - int(int(data["Nodes_X"]) / 2)) * data["Nodes_X"]
        self.pos_x = x  
        
    def generateVolume1D(data,horizon_factor):
        volumes = np.empty(data["Nodes_X"])
        volumes.fill(float(data["Length_X"])*float(data["Surface_X"]))
       
        for i in range(0, horizon_factor):
            volume_boundary += volumes[i]

        for i in range(len(self.volumes) - horizon_factor, len(self.volumes) ):
            volume_boundary += volumes[i]
       
    #def generateGrid2D():
        
    #def generateGrid3D():
    
    
    #def read_nodes():
        
        
    #def read_bc():
