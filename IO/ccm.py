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
        self.y = problem.y
        self.force_int = problem.force_int
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
        
        

    