#-*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca

from ..util import linalgebra
import numpy as np
from ..util import functions


## Abstract class of the problem classes, which contains common methods 
class Problem():
    
    # @param deck The input deck
    def compute_volume_correction(self,deck):
        ## Volume correction factor for each node
        self.volume_correction = np.ones( ( deck.num_nodes, self.neighbors.max_neighbors), dtype=np.float64 )
        for i in range(0, deck.num_nodes):
            index_x_family = self.neighbors.get_index_x_family(i)
            n = 0
            for p in index_x_family:
                X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                r = deck.delta_X / 2.0
                if linalgebra.norm(X) > self.neighbors.horizon - r:
                    self.volume_correction[i,n] = (self.neighbors.horizon + r - linalgebra.norm(X)) / (deck.delta_X)
                else:
                    pass
                n += 1

    ## Compute the weighted volume for each node
    # @param deck The input deck
    def compute_weighted_volume(self, deck):
        ## Weighted volume for each node
        self.weighted_volume = np.zeros((deck.num_nodes),dtype=np.float64)
        for i in range(0, deck.num_nodes):
            index_x_family = self.neighbors.get_index_x_family(i)
            n = 0
            for p in index_x_family:
                X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                self.weighted_volume[i] += functions.w(deck, X, deck.influence_function) * (linalgebra.norm(X))**2 * self.volume_correction[i,n] * deck.geometry.volumes[p]
                n += 1    