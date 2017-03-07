# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 16:16:07 2015

@author: ilyass.tabiai@polymtl.ca
@author: rolland.delorme@polymtl.ca
@author: patrick.diehl@polymtl.ca
"""

import numpy as np

class elastic_material():
    @profile
    def __init__(self, PD_deck, PD_problem, y):
        self.len_x = PD_deck.num_nodes_x
        self.E_Modulus = PD_deck.e_modulus
        self.compute_ext_state(PD_deck, PD_problem, y)
        self.compute_T(PD_deck, PD_problem, y)
        self.compute_Ts(PD_deck, PD_problem)

    # Computes the deformations using the current positions and initial
    # positions of each node
    # Reminder: y[xi, t] = x[xi] + u[xi, t]
    @profile
    def compute_ext_state(self, PD_deck, PD_problem, y):
        # Initialization for e
        e = np.zeros((self.len_x, self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family(x_i)
            for x_p in index_x_family:
                e[x_i, x_p] = np.absolute(y[x_p] - y[x_i]) - np.absolute(PD_deck.geometry.pos_x[x_p] - PD_deck.geometry.pos_x[x_i])
        self.e = e
    @profile
    def compute_T(self, PD_deck, PD_problem, y):
        w = PD_deck.influence_function
        M = PD_problem.compute_m(y)
        tscal = np.zeros((self.len_x, self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family(x_i)
            for x_p in index_x_family:
                tscal[x_i, x_p] = (w / PD_problem.weighted_function(PD_deck,PD_deck.geometry.pos_x,x_i)) * self.E_Modulus * self.e[x_i,x_p]
        self.tscal = tscal

        T = np.zeros((self.len_x, self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family(x_i)
            for x_p in index_x_family:
                T[x_i, x_p] = tscal[x_i, x_p] * M[x_i, x_p]
        self.T = T
    @profile
    def compute_Ts(self, PD_deck, PD_problem):
        Ts = np.zeros((self.len_x))
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family(x_i)
            for x_p in index_x_family:
                Ts[x_i] = Ts[x_i] + self.T[x_i, x_p] - self.T[x_p, x_i]
            Ts[x_i] = Ts[x_i] * PD_deck.geometry.volumes[x_i]
        self.Ts = Ts
