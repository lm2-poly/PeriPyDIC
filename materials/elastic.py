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
        ## Scalar influence function
        self.w = deck.influence_function
        if deck.dim == 1:
            ## Young modulus of the material
            self.Young_Modulus = deck.young_modulus
            ##  Deformed direction vector state in x direction
            self.M_x = problem.compute_m(y,deck.dim,deck.num_nodes)
        if deck.dim == 2:
            ## Bulk modulus of the material
            self.K = deck.bulk_modulus
            ## Shear modulus of the material
            self.Mu = deck.shear_modulus
            ## Poisson ratio of the material
            self.Nu = (3 * self.K - 2 * self.Mu) / (2 * (3 * self.K + self.Mu))
            ##  Deformed direction vector state in y direction
            self.M_x , self.M_y = problem.compute_m(y,deck.dim,deck.num_nodes)
        if deck.dim == 3:
            ## Bulk modulus of the material
            self.K = deck.bulk_modulus
            ## Shear modulus of the material
            self.Mu = deck.shear_modulus
            ##  Deformed direction vector state in z direction
            self.M_x , self.M_y , self.M_z = problem.compute_m(y,deck.dim,deck.num_nodes)
            
        self.compute_ext_state(deck, problem, y)
        self.compute_T(deck, problem, y)
        self.compute_Ts(deck, problem)

    # Computes the dilatation
    def pd_dilatation(self, deck, problem, x, e, x_i):
        index_x_family = problem.neighbors.get_index_x_family(x_i)
        result = 0
        dilatation = 0
        for x_p in index_x_family:
            
            if deck.dim == 1:
                actual = np.absolute(x[x_p] - x[x_i]) * e[x_i,x_p]
                result += deck.influence_function * actual * deck.geometry.volumes[x_p]
                dilatation = (1 / problem.weighted_function(deck,deck.geometry.nodes,x_i)) * result 
            
            if deck.dim == 2:
                actual = np.sqrt(np.power(x[x_p][0] - x[x_i][0],2) + np.power(x[x_p][1] - x[x_i][1],2)) * e[x_i, x_p]
                result += deck.influence_function * actual * deck.geometry.volumes[x_p]
                dilatation = (2 / problem.weighted_function(deck,deck.geometry.nodes,x_i)) * ((2 * self.Nu - 1) / (self.Nu - 1)) * result
            
            if deck.dim == 3:
                actual = np.sqrt(np.power(x[x_p][0] - x[x_i][0],2) + np.power(x[x_p][1] - x[x_i][1],2) + np.power(x[x_p][2] - x[x_i][2],2)) * e[x_i,x_p]
                result += deck.influence_function * actual * deck.geometry.volumes[x_p]
                dilatation = (3 / problem.weighted_function(deck,deck.geometry.nodes,x_i)) * result 
        return dilatation
  
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
                    actual = np.absolute(y[x_p] - y[x_i])
                    initial = np.absolute(deck.geometry.nodes[x_p] - deck.geometry.nodes[x_i])
                    e[x_i, x_p] = actual - initial
                    self.e = e
                
                if deck.dim == 2:
                    e_s = np.zeros((self.len_x, self.len_x))
                    e_d = np.zeros((self.len_x, self.len_x))
                    actual = np.sqrt(np.power(y[x_p] - y[x_i],2) + np.power(y[self.len_x + x_p] - y[self.len_x + x_i],2))
                    initial = np.sqrt(np.power(deck.geometry.nodes[x_p][0] - deck.geometry.nodes[x_i][0],2)+np.power(deck.geometry.nodes[x_p][1] - deck.geometry.nodes[x_i][1],2))
                    e[x_i, x_p] = actual  - initial
                    e_s[x_i, x_p] = self.pd_dilatation(deck, problem, deck.geometry.nodes, e, x_i) * initial / 3
                    e_d[x_i, x_p] = e[x_i, x_p] - e_d[x_i, x_p]
                    ## Scalar extension state
                    self.e = e 
                    self.e_s = e_s
                    self.e_d = e_d
                    
                if deck.dim == 3:
                    e_s = np.zeros((self.len_x, self.len_x))
                    e_d = np.zeros((self.len_x, self.len_x))
                    actual = np.sqrt(np.power(y[x_p] - y[x_i],2) + np.power(y[self.len_x + x_p] - y[self.len_x + x_i],2)+ np.power(y[2*self.len_x + x_p] - y[2*self.len_x + x_i],2))
                    initial = np.sqrt(np.power(deck.geometry.nodes[x_p][0] - deck.geometry.nodes[x_i][0],2)+np.power(deck.geometry.nodes[x_p][1] - deck.geometry.nodes[x_i][1],2)+np.power(deck.geometry.nodes[x_p][2] - deck.geometry.nodes[x_i][2],2))
                    e[x_i, x_p] = actual  - initial
                    e_s[x_i, x_p] = self.pd_dilatation(deck, problem, deck.geometry.nodes, e, x_i) * initial / 3
                    e_d[x_i, x_p] = e[x_i, x_p] - e_d[x_i, x_p]
                    ## Scalar extension state
                    self.e = e 
                    self.e_s = e_s
                    self.e_d = e_d
    
    ## Function to compute the vector force state
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual postions
    def compute_T(self, deck, problem, y):
        tscal = np.zeros((self.len_x, self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                
                if deck.dim == 1:
                    alpha = self.Young_Modulus / problem.weighted_function(deck,deck.geometry.nodes,x_i)
                    tscal[x_i, x_p] = alpha * self.w * self.e[x_i,x_p]
                
                if deck.dim == 2:
                    tscal_s = np.zeros((self.len_x, self.len_x))
                    tscal_d = np.zeros((self.len_x, self.len_x))
                    alpha_s = (9 / problem.weighted_function(deck,deck.geometry.nodes,x_i)) * self.K
                    alpha_d = (8 / problem.weighted_function(deck,deck.geometry.nodes,x_i)) * self.Mu
                    tscal_s[x_i, x_p] = alpha_s * self.w * self.e_s[x_i,x_p]
                    tscal_d[x_i, x_p] = alpha_d * self.w * self.e_d[x_i,x_p]
                    tscal[x_i, x_p] = tscal_s[x_i, x_p] + tscal_d[x_i, x_p]
                
                if deck.dim == 3:
                    tscal_s = np.zeros((self.len_x, self.len_x))
                    tscal_d = np.zeros((self.len_x, self.len_x))
                    alpha_s = (9 / problem.weighted_function(deck,deck.geometry.nodes,x_i)) * self.K
                    alpha_d = (15 / problem.weighted_function(deck,deck.geometry.nodes,x_i)) * self.K
                    tscal_s[x_i, x_p] = alpha_s * self.w * self.e_s[x_i,x_p]
                    tscal_d[x_i, x_p] = alpha_d * self.w * self.e_d[x_i,x_p]
                    tscal[x_i, x_p] = tscal_s[x_i, x_p] + tscal_d[x_i, x_p]
                    
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
