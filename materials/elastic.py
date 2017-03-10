# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import numpy as np

## Class to compute, as a vector state, the global internal volumic force of an elastic material using its material properties
class Elastic_material():
  
    def __init__(self, PD_deck, PD_problem, y):
        ## Amount of nodes in x direction
        self.len_x = PD_deck.num_nodes_x
        ## Young's modulus of the material
        self.E_Modulus = PD_deck.e_modulus
        ## Scalar influence function
        self.w = PD_deck.influence_function
        ##  Deformed direction vector state
        self.M = PD_problem.compute_m(y)
        self.compute_ext_state(PD_deck, PD_problem, y)
        self.compute_T(PD_deck, PD_problem, y)
        self.compute_Ts(PD_deck, PD_problem)

    # Computes the deformations using the current positions and initial
    # positions of each node
    # Reminder: y[xi, t] = x[xi] + u[xi, t]
  
    ## Function to compute the scalar extension state
    # @param PD_deck blabla...
    def compute_ext_state(self, PD_deck, PD_problem, y):
        # Initialization for e
        e = np.zeros((self.len_x, self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family(x_i)
            for x_p in index_x_family:
                e[x_i, x_p] = np.absolute(y[x_p] - y[x_i]) - np.absolute(PD_deck.geometry.pos_x[x_p] - PD_deck.geometry.pos_x[x_i])
        ## Scalar extension state
        self.e = e
    
    ## Function to compute the vector force state
    def compute_T(self, PD_deck, PD_problem, y):
        tscal = np.zeros((self.len_x, self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family(x_i)
            for x_p in index_x_family:
                tscal[x_i, x_p] = (self.w / PD_problem.weighted_function(PD_deck,PD_deck.geometry.pos_x,x_i)) * self.E_Modulus * self.e[x_i,x_p]
        ## Scalar force state
        self.tscal = tscal

        T = np.zeros((self.len_x, self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family(x_i)
            for x_p in index_x_family:
                T[x_i, x_p] = tscal[x_i, x_p] * self.M[x_i, x_p]
        ##  Vector force state
        self.T = T
    
    ## Function to compute, as a vector state, the global internal volumic force within the equation of motion
    def compute_Ts(self, PD_deck, PD_problem):
        Ts = np.zeros((self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family(x_i)
            for x_p in index_x_family:
                Ts[x_i] = Ts[x_i] + self.T[x_i, x_p] - self.T[x_p, x_i]
            Ts[x_i] = Ts[x_i] * PD_deck.geometry.volumes[x_i]
        ## Sum of (T[xi] - T[xp])*Vol[xi] equivalent to the global internal volumic force
        self.Ts = Ts
