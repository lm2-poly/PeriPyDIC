# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import numpy as np
from scipy import linalg
import sys
np.set_printoptions(threshold='nan')

## Class to compute, as a vector state, the global internal volumic force of an elastic material using its material properties
class CCM_calcul():

    ## Constructor
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual postions
    def __init__(self, deck, problem):
        ## Scalar influence function
        self.x = deck.geometry.nodes
        self.y = problem.y
        self.force_int = problem.force_int
        self.ext = problem.ext
        self.num_nodes = deck.num_nodes
        self.dim = deck.dim 
        self.time_steps = deck.time_steps
        
        self.compute_u_displacement()

    # Compute the displacement for each Node
    def compute_u_displacement(self):
        self.u = np.zeros((self.num_nodes, self.dim, self.time_steps),dtype=np.float64)
        for t_n in range(1, self.time_steps):
            for i in range(0, self.num_nodes):
                self.u[i,:,t_n] = self.y[i,:,t_n] - self.y[i,:,t_n-1]

            
    # Compute the the global internal force density vector
    def X_vector_state(self, i, p):
        index_x_family = problem.neighbors.get_index_x_family(i)
        if p in index_x_family:
            X = self.x[p,:] - self.x[i,:]
        else:
            X = np.zeros((1, self.dim),dtype=np.float64)
            