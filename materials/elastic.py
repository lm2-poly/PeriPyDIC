# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import numpy as np
from scipy import linalg
import sys
np.set_printoptions(threshold='nan')

## Class to compute, as a vector state, the global internal volumic force of an elastic material using its material properties
class Elastic_material():

    ## Constructor
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual postions
    def __init__(self, deck, problem, y):
        ## Scalar influence function
        self.w = deck.influence_function

        ## Weighted volume
        self.Weighted_Volume = problem.weighted_volume
        
        if deck.dim == 1:
            ## Young modulus of the material
            self.Young_Modulus = deck.young_modulus
        if deck.dim == 2:
            ## Bulk modulus of the material
            self.K = deck.bulk_modulus
            ## Shear modulus of the material
            self.Mu = deck.shear_modulus
            ## Poisson ratio of the material
            self.Nu = (3. * self.K - 2. * self.Mu) / (2. * (3. * self.K + self.Mu))

        if deck.dim == 3:
            ## Bulk modulus of the material
            self.K = deck.bulk_modulus
            ## Shear modulus of the material
            self.Mu = deck.shear_modulus

        self.compute_dilatation(deck, problem, y)
        self.compute_f_int(deck, problem, y)

    # Compute the dilatation for each Node
    def compute_dilatation(self, deck, problem, y):
        self.dilatation = np.zeros((deck.num_nodes),dtype=np.float64)
        self.e = np.zeros((deck.num_nodes, deck.num_nodes),dtype=np.float64)
        for i in range(0, deck.num_nodes):
            index_x_family = problem.neighbors.get_index_x_family(i)
            for p in index_x_family:
                    Y = (y[p,:]) - y[i,:]
                    X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                    self.e[i,p] = linalg.norm(Y) - linalg.norm(X)
                    
                    if deck.dim == 1:
                        self.dilatation[i] += (1. / self.Weighted_Volume[i]) * self.w * linalg.norm(X) * self.e[i,p] * deck.geometry.volumes[p]
        
                    if deck.dim == 2:
                        self.dilatation[i] += (2. / self.Weighted_Volume[i]) * ((2. * self.Nu - 1.) / (self.Nu - 1.)) * self.w * linalg.norm(X) * self.e[i,p] * deck.geometry.volumes[p]
        
                    if deck.dim == 3:
                        self.dilatation[i] += (3. / self.Weighted_Volume[i]) * self.w * linalg.norm(X) * self.e[i,p] * deck.geometry.volumes[p]

    # Compute the the global internal force density vector
    def compute_f_int(self, deck, problem, y):
        self.f_int = np.zeros((deck.num_nodes, deck.dim),dtype=np.float64)
        for i in range(0, deck.num_nodes):
            index_x_family = problem.neighbors.get_index_x_family(i)
            for p in index_x_family:
                Y = y[p,:] - y[i,:]
                X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                
                # Compute the direction vector between Node_p and Node_i
                M = Y / linalg.norm(Y)

                if deck.dim == 1:
                    # PD material parameter
                    alpha = self.Young_Modulus / self.Weighted_Volume[i]
                    # Scalar force state
                    self.t = alpha * self.w * self.e[i,p]

                if deck.dim == 2:
                    # PD material parameter
                    alpha_s = (9. / self.Weighted_Volume[i]) * self.K
                    alpha_d = (8. / self.Weighted_Volume[i]) * self.Mu
                    # Scalar force state
                    e_s = self.dilatation[i] * linalg.norm(X) / 3.
                    e_d = self.e[i, p] - e_s
                    t_s = alpha_s * self.w * e_s
                    t_d = alpha_d * self.w * e_d
                    self.t = t_s + t_d

                if deck.dim == 3:
                    # PD material parameter
                    alpha_s = (9. / self.Weighted_Volume[i]) * self.K
                    alpha_d = (15. / self.Weighted_Volume[i]) * self.K
                    # Scalar force state
                    e_s = self.dilatation[i] * linalg.norm(X) / 3.
                    e_d = self.e[i, p] - e_s
                    t_s = alpha_s * self.w * e_s
                    t_d = alpha_d * self.w * e_d
                    self.t = t_s + t_d
                
                self.f_int[i,:] += self.t * M * deck.geometry.volumes[p]
                self.f_int[p,:] += -self.t * M * deck.geometry.volumes[i]
                
