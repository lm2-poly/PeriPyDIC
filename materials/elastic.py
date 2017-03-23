# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import numpy as np
from numpy import linalg
import sys
np.set_printoptions(threshold='nan')

## Class to compute, as a vector state, the global internal volumic force of an elastic material using its material properties
class Elastic_material():

    ## Constructor
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual postions
    def __init__(self, deck, problem, y, pertub=0.0):
        ## Scalar influence function
        self.w = deck.influence_function

        ## Weighted volume
        self.Weighted_Volume = problem.weighted_volume
        
        self.pertub = pertub

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

        self.compute_ext_state(deck, problem, y)
        self.compute_tscal(deck, problem)
        self.compute_T(deck, problem, y)
        self.compute_F(deck, problem)

    # Computes the dilatation at Node_i
    def pd_dilatation(self, deck, problem, e, i):
        dilatation = 0.
        index_x_family = problem.neighbors.get_index_x_family(i)
        for p in index_x_family:
            if deck.dim == 1:
                X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                dilatation += (1. / self.Weighted_Volume[i]) * self.w * np.linalg.norm(X) * e[i,p] * deck.geometry.volumes[p]

            if deck.dim == 2:
                X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                dilatation += (2. / self.Weighted_Volume[i]) * ((2. * self.Nu - 1.) / (self.Nu - 1.)) * self.w * np.linalg.norm(X) * e[i,p] * deck.geometry.volumes[p]

            if deck.dim == 3:
                X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                dilatation += (3. / self.Weighted_Volume[i]) * self.w * np.linalg.norm(X) * e[i,p] * deck.geometry.volumes[p]
        return dilatation

    # Computes the direction vector between Node_p and Node_i
    def compute_dir_vector(self, deck, y, i, p):
        Y = (y[p,:]+self.pertub) - y[i,:]
        M = Y / np.linalg.norm(Y)
        return M

    ## Compute the scalar extension state between Node_p and Node_i for each node
    def compute_ext_state(self, deck, problem, y):
        # Initialization for e
        self.e = np.zeros((deck.num_nodes, deck.num_nodes),dtype=np.float64)
        for i in range(0, deck.num_nodes):
            index_x_family = problem.neighbors.get_index_x_family(i)
            for p in index_x_family:
                if deck.dim >=1:
                    Y = (y[p,:]+self.pertub) - y[i,:]
                    X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                    self.e[i,p] = np.linalg.norm(Y) - np.linalg.norm(X)

                if deck.dim >= 2:
                    self.e_s = np.zeros((deck.num_nodes, deck.num_nodes),dtype=np.float64)
                    self.e_d = np.zeros((deck.num_nodes, deck.num_nodes),dtype=np.float64)
                    self.e_s[i, p] = self.pd_dilatation(deck, problem, self.e, i) * np.linalg.norm(X) / 3.
                    self.e_d[i, p] = self.e[i, p] - self.e_s[i, p]

    ## Compute the scalar force state between between Node_p and Node_i for each node
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual postions
    def compute_tscal(self, deck, problem):
        self.tscal = np.zeros((deck.num_nodes, deck.num_nodes),dtype=np.float64)
        for i in range(0, deck.num_nodes):
            index_x_family = problem.neighbors.get_index_x_family(i)
            for p in index_x_family:

                if deck.dim == 1:
                    # PD material parameter
                    alpha = self.Young_Modulus / self.Weighted_Volume[i]
                    # Scalar force state
                    self.tscal[i,p] = alpha * self.w * self.e[i,p]

                if deck.dim == 2:
                    # PD material parameter
                    alpha_s = (9. / self.Weighted_Volume[i]) * self.K
                    alpha_d = (8. / self.Weighted_Volume[i]) * self.Mu
                    # Scalar force state
                    self.tscal_s = np.zeros((deck.num_nodes, deck.num_nodes),dtype=np.float64)
                    self.tscal_d = np.zeros((deck.num_nodes, deck.num_nodes),dtype=np.float64)
                    self.tscal_s[i,p] = alpha_s * self.w * self.e_s[i,p]
                    self.tscal_d[i,p] = alpha_d * self.w * self.e_d[i,p]
                    self.tscal[i,p] = self.tscal_s[i, p] + self.tscal_d[i,p]

                if deck.dim == 3:
                    # PD material parameter
                    alpha_s = (9. / self.Weighted_Volume[i]) * self.K
                    alpha_d = (15. / self.Weighted_Volume[i]) * self.K
                    # Scalar force state
                    self.tscal_s = np.zeros((deck.num_nodes, deck.num_nodes),dtype=np.float64)
                    self.tscal_d = np.zeros((deck.num_nodes, deck.num_nodes),dtype=np.float64)
                    self.tscal_s[i,p] = alpha_s * self.w * self.e_s[i,p]
                    self.tscal_d[i,p] = alpha_d * self.w * self.e_d[i,p]
                    self.tscal[i,p] = self.tscal_s[i,p] + self.tscal_d[i,p]

    ## Compute the vector force state between between Node_p and Node_i for each node
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual positions
    def compute_T(self, deck, problem, y):
        self.T = np.zeros((deck.num_nodes, deck.dim, deck.num_nodes),dtype=np.float64)
        for i in range(0, deck.num_nodes):
            index_x_family = problem.neighbors.get_index_x_family(i)
            for p in index_x_family:
                self.T[i,:,p] = self.tscal[i,p] * self.compute_dir_vector( deck, y, i, p)

    ## Compute the global vector force state for the equation of motion
    # @param deck The input deck
    # @param problem The related peridynamic problem
    def compute_F(self, deck, problem):
        self.F = np.zeros((deck.num_nodes, deck.dim),dtype=np.float64)
        for i in range(0, deck.num_nodes):
            index_x_family = problem.neighbors.get_index_x_family(i)
            for p in index_x_family:
                self.F[i,:] += (self.T[i,:,p] - self.T[p,:,i]) * deck.geometry.volumes[p]
        #print self.F
