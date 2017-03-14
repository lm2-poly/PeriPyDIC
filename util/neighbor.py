#-*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import numpy as np
import scipy.spatial
import sys

## Class for handling the neighborhood
class NeighborSearch():
    
    ## Constructor
    # @param The input deck
    def __init__(self,deck):
        ## Saftey factor for the search of the neighborhood
        # @param deck The input deck
        self.safety_factor = 1.001
        ## Horizon of the neighborhood
        self.horizon = deck.horizon_factor_m_value*deck.delta_x*self.safety_factor
        self.findNeighbors(deck)
        
    ## Returns the family of node x_i
    # @param x_i Id of the node
    # @return The ids of the neighbors of node x_i
    def get_index_x_family(self, x_i):
        return self.family[x_i]

    ## Generates adjacency matrix
    #@param The input deck
    def findNeighbors(self,deck):
        tree = scipy.spatial.KDTree(deck.geometry.nodes)    
        d , ids = families = tree.query(deck.geometry.nodes, k=20, eps=0.0, p=2, distance_upper_bound=self.horizon) 
        self.family = []
        for i in range(0,len(ids)):
           self.family.append(np.array(filter(lambda x: x != i, np.array(filter(lambda x: x < deck.num_nodes, ids[i])))))
        del d , ids
        
