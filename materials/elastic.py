# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import numpy as np

## Class to compute, as a vector state, the global internal volumic force of an elastic material using its material properties
class Elastic_material():
  
    ## Constructor
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual postions
    def __init__(self, deck, problem, y):
        ## Amount of nodes in x direction
        self.len_x = deck.num_nodes
        ## Young's modulus of the material
        self.E_Modulus = deck.e_modulus
        ## Scalar influence function
        self.w = deck.influence_function
        ##  Deformed direction vector state
        self.M = problem.compute_m(y)
        self.compute_ext_state(deck, problem, y)
        self.compute_T(deck, problem, y)
        self.compute_Ts(deck, problem)

    # Computes the deformations using the current positions and initial
    # positions of each node
    # Reminder: y[xi, t] = x[xi] + u[xi, t]
  
    ## Function to compute the scalar extension state
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual postions
    def compute_ext_state(self, deck, problem, y):
        # Initialization for e
        e = np.zeros((self.len_x, self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                e[x_i, x_p] = np.absolute(y[x_p] - y[x_i]) - np.absolute(deck.geometry.nodes[x_p] - deck.geometry.nodes[x_i])
        ## Scalar extension state
        self.e = e
    
    ## Function to compute the vector force state
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual postions
    def compute_T(self, deck, problem, y):
        tscal = np.zeros((self.len_x, self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                tscal[x_i, x_p] = (self.w / problem.weighted_function(deck,deck.geometry.nodes,x_i)) * self.E_Modulus * self.e[x_i,x_p]
        ## Scalar force state
        self.tscal = tscal

        T = np.zeros((self.len_x, self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                T[x_i, x_p] = tscal[x_i, x_p] * self.M[x_i, x_p]
        ##  Vector force state
        self.T = T
    
    ## Function to compute, as a vector state, the global internal volumic force within the equation of motion
    # @param deck The input deck
    # @param problem The related peridynamic problem
    def compute_Ts(self, deck, problem):
        Ts = np.zeros((self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                Ts[x_i] = Ts[x_i] + self.T[x_i, x_p] - self.T[x_p, x_i]
            Ts[x_i] = Ts[x_i] * deck.geometry.volumes[x_i]
        ## Sum of (T[xi] - T[xp])*Vol[xi] equivalent to the global internal volumic force
        self.Ts = Ts
