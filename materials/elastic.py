# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import numpy as np
#import sys
np.set_printoptions(threshold='nan')

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
        if deck.dim == 1:
            ##  Deformed direction vector state in x direction
            self.M_x = problem.compute_m(y,deck.dim,deck.num_nodes)
        if deck.dim == 2:
            ##  Deformed direction vector state in y direction
            self.M_x , self.M_y = problem.compute_m(y,deck.dim,deck.num_nodes)
        if deck.dim == 3:
            ##  Deformed direction vector state in z direction
            self.M_x , self.M_y , self.M_z = problem.compute_m(y,deck.dim,deck.num_nodes)
            
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
                if deck.dim == 1:
                    e[x_i, x_p] = np.absolute(y[x_p] - y[x_i]) - np.absolute(deck.geometry.nodes[x_p] - deck.geometry.nodes[x_i])
                if deck.dim == 2:
                    actual = np.sqrt(np.power(y[x_p] - y[x_i],2) + np.power(y[self.len_x + x_p] - y[self.len_x + x_i],2))
                    initial = np.sqrt(np.power(deck.geometry.nodes[x_p][0] - deck.geometry.nodes[x_i][0],2)+np.power(deck.geometry.nodes[x_p][1] - deck.geometry.nodes[x_i][1],2))
                    e[x_i, x_p] = actual  - initial 
                if deck.dim == 3:
                    actual = np.sqrt(np.power(y[x_p] - y[x_i],2) + np.power(y[self.len_x + x_p] - y[self.len_x + x_i],2)+ np.power(y[2*self.len_x + x_p] - y[2*self.len_x + x_i],2))
                    initial = np.sqrt(np.power(deck.geometry.nodes[x_p][0] - deck.geometry.nodes[x_i][0],2)+np.power(deck.geometry.nodes[x_p][1] - deck.geometry.nodes[x_i][1],2)+np.power(deck.geometry.nodes[x_p][2] - deck.geometry.nodes[x_i][2],2))
                    e[x_i, x_p] = actual  - initial 
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
        ## Vector state in x direction
        self.T_x = np.zeros((deck.dim * self.len_x, deck.dim * self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                if deck.dim >= 1:
                    self.T_x[x_i, x_p] = tscal[x_i, x_p] * self.M_x[x_i, x_p]
                if deck.dim >=2:
                    ## Vector state in y direction
                    self.T_y = np.zeros((deck.dim * self.len_x, deck.dim * self.len_x))
                    self.T_y[x_i, x_p] = tscal[x_i, x_p] * self.M_y[x_i, x_p]
                if deck.dim >= 3:
                    ## Vector state in z direction
                    self.T_z = np.zeros((deck.dim * self.len_x, deck.dim * self.len_x))
                    self.T_z[x_i, x_p] = tscal[x_i, x_p] * self.M_z[x_i, x_p]
        ##  Vector force state
       
    
    
    ## Function to compute, as a vector state, the global internal volumetric force within the equation of motion
    # @param deck The input deck
    # @param problem The related peridynamic problem
    def compute_Ts(self, deck, problem):
        Ts = np.zeros((deck.dim *self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                if deck.dim >= 1:
                    Ts[x_i] = Ts[x_i] + self.T_x[x_i, x_p] - self.T_x[x_p, x_i]
                if deck.dim >= 2:
                    Ts[self.len_x + x_i] = Ts[self.len_x + x_i] + self.T_y[x_i, x_p] - self.T_y[x_p, x_i] 
                if deck.dim >= 3:
                    Ts[2*self.len_x + x_i] = Ts[2*self.len_x + x_i] + self.T_z[x_i, x_p] - self.T_z[x_p, x_i]        
            
            
            if deck.dim >= 1:
                Ts[x_i] = Ts[x_i] * deck.geometry.volumes[x_i]
            if deck.dim >= 2:
                Ts[self.len_x + x_i] = Ts[self.len_x + x_i] * deck.geometry.volumes[x_i]
            if deck.dim >= 3:
                Ts[2*self.len_x + x_i] = Ts[2*self.len_x + x_i] * deck.geometry.volumes[x_i]
        ## Sum of (T[xi] - T[xp])*Vol[xi] equivalent to the global internal volumic force
        self.Ts = Ts
