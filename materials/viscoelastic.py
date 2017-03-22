# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca
import numpy as np

## Class to compute, as a vector state, the global internal volumic force of an viscoelastic material using its material properties
class Viscoelastic_material():

    ## Constructor
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual postions
    # @param t_n The time step
    def __init__(self, deck, problem, y, t_n):
        ## Amount of nodes in x direction
        self.len_x = deck.num_nodes
        ## Relaxation modulus of the material
        self.Relax_Modulus = deck.relax_modulus
        ## Relaxation time of the material
        self.Relax_Time = deck.relax_time
        ## Scalar influence function
        self.w = deck.influence_function
        ##  Deformed direction vector state
        self.M = problem.compute_m(y)
        self.compute_ext_state(deck, problem, y)
        self.compute_ext_state_visco(deck, problem, y, t_n )
        self.compute_t_visco(deck, problem, y)
        self.compute_T(deck, problem, y)
        self.compute_Ts(deck, problem)

    #Computes the deformations using the current positions and initial
    #positions of each node
    #Reminder: y[xi, t] = x[xi] + u[xi, t]

    ## Function to compute the scalar extension state
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual postions
    def compute_ext_state(self, deck, problem, y):
        e = np.zeros( (self.len_x, self.len_x) )
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                e[x_i, x_p] = np.absolute(y[x_p] - y[x_i]) - np.absolute(deck.geometry.nodes[x_p] - deck.geometry.nodes[x_i])
        ## Scalar extension state
        self.e = e

    ## Function to compute the viscoelastic part of the scalar extension state
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual postions
    # @param t_n The time step
    def compute_ext_state_visco(self, deck, problem, y, t_n):

        e_visco = np.zeros( (self.len_x, self.len_x, len(self.Relax_Time)) )
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                for k in range(1, len(self.Relax_Time)):
                    temp_exp = np.exp((- deck.delta_t)/(self.Relax_Time[k]))
                    delta_e = self.e[x_i, x_p] - problem.ext[x_i, x_p, t_n-1]
                    beta = 1-(self.Relax_Time[k]*(1-temp_exp)) / deck.delta_t
                    e_visco[x_i, x_p, k] = problem.ext[x_i, x_p, t_n-1] *(1-temp_exp) + problem.ext_visco[x_i, x_p, k, t_n-1]*temp_exp + beta*delta_e
        ## Viscous part of the scalar extension state
        self.e_visco = e_visco

    ## Function to compute the viscous part of the scalar force state
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual postions
    def compute_t_visco(self, deck, problem, y):
        t_visco = np.zeros( (self.len_x, self.len_x ) )
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                for k in range(1, len(self.Relax_Time)):
                    t_visco[x_i, x_p] = t_visco[x_i, x_p] + (self.w / problem.weighted_function(deck, deck.geometry.nodes, x_i))*self.Relax_Modulus[k]*(self.e[x_i, x_p] - self.e_visco[x_i, x_p, k])
        ## Viscous part of the scalar state
        self.t_visco = t_visco

    ## Function to compute the vector force state
    # @param deck The input deck
    # @param problem The related peridynamic problem
    # @param y The actual postions
    def compute_T(self, deck, problem, y):
        tscal = np.zeros( (self.len_x, self.len_x ) )
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                tscal[x_i, x_p] = (self.w / problem.weighted_function(deck, deck.geometry.nodes, x_i))*self.Relax_Modulus[0]*self.e[x_i, x_p] + self.t_visco[x_i, x_p]
        ## Scalar force state
        self.tscal = tscal

        T = np.zeros( (self.len_x, self.len_x ) )
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                T[x_i, x_p] = tscal[x_i, x_p] * self.M[x_i, x_p]
        ## Vector force state
        self.T = T

    ## Function to compute, as a vector state, the global internal volumic force within the equation of motion
    # @param deck The input deck
    # @param problem The related peridynamic problem
    def compute_Ts(self, deck, problem):
        Ts = np.zeros( (self.len_x ) )
        for x_i in range(0, self.len_x):
            index_x_family = problem.neighbors.get_index_x_family(x_i)
            for x_p in index_x_family:
                Ts[x_i] = Ts[x_i] + self.T[x_i, x_p] - self.T[x_p, x_i]
            Ts[x_i] = Ts[x_i] * deck.geometry.volumes[x_i]
        ## Sum of the (T[xi] - T[xp])*Vol[xi] equivalent to the global internal volumic force
        self.Ts = Ts
