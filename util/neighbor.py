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
        if deck.geometry.dim == 1:
            self.generate_neighborhood_matrix(deck)
        else:
            self.findNeighbors(deck)
        
    ## Returns the family of node x_i
    # @param x_i Id of the node
    # @return The ids of the neighbors of node x_i
    def get_index_x_family(self, x_i):
        return (np.where(self.family[x_i] == 1))[0]

    ## Generates adjacency matrix
    #@param The input deck
    def generate_neighborhood_matrix(self, deck):
        ## Adjacency matrix 
        self.family = np.zeros((deck.num_nodes, deck.num_nodes))
        for x_i in range(0, deck.num_nodes):
            for x_p in range(0, deck.num_nodes):
                if x_p == x_i:
                    pass
                elif np.absolute(deck.geometry.nodes[x_i] - deck.geometry.nodes[x_p]) <= self.horizon:
                    self.family[x_i][x_p] = 1
                else:
                    pass
        #print self.family
        
    def findNeighbors(self,deck):
        tree = scipy.spatial.cKDTree(deck.geometry.nodes)    
        print deck.geometry.nodes
        deq, ieq = tree.query(deck.geometry.nodes, k=10, eps=0.0, p=2, distance_upper_bound=self.horizon) 
        print deq
        sys.exit()
        
        
