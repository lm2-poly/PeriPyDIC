#-*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import numpy as np
import scipy.spatial
import sys 

sys.setrecursionlimit(100000)

## Class for handling the neighborhood
class NeighborSearch():

    ## Constructor
    # @param deck The input deck
    def __init__(self,deck):
        ## Safety factor for the search of the neighborhood
        self.safety_factor = deck.safety_factor
        ## Maximal amount of neighbors
        self.max_neighbors = 0
        ## Horizon of the neighborhood
        self.horizon = deck.horizon_factor_m_value * deck.delta_X * self.safety_factor
        self.findNeighbors(deck)
        
        if deck.critical_stretch > 0:
            self.damage = np.zeros((deck.num_nodes,deck.num_nodes))

    ## Returns the family of node "i"
    # @param i Id of the node
    # @return The ids of the neighbors of node "i"
    def get_index_x_family(self, i):
        return self.family[i]

    ## Generates adjacency lists
    # @param deck The input deck
    def findNeighbors(self,deck):
        tree = scipy.spatial.KDTree(deck.geometry.nodes)
        d , ids =  tree.query(deck.geometry.nodes, k=150, eps=0.0, p=2, distance_upper_bound=self.horizon)
        ## Adjaceny list
        self.family = []
        for i in range(0,len(ids)):
            self.family.append(np.array(filter(lambda x: x != i, np.array(filter(lambda x: x < deck.num_nodes, ids[i])))))
            tmp = len(self.family[i])
            if tmp > self.max_neighbors:
                self.max_neighbors = tmp
        del d , ids
