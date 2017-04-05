# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 10:52:04 2016

@author: ilyass
"""
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 16:16:07 2015

@author: ilyass.tabiai@gmail.com
@author: rolland.delorme@gmail.com
@author: diehl@ins.uni-bonn.de
"""

from problem import PD_problem
import numpy as np

class elastic_material_dic():

    def __init__(self, PD_deck, PD_problem, t):
        self.Exp_Num_Nodes = len(PD_problem.exp_displacement[0])
        self.Modulus = PD_deck.get_elastic_material_properties()

        self.exp_Ts = np.zeros(int(self.Exp_Num_Nodes))
        self.exp_W = np.zeros(int(self.Exp_Num_Nodes))

        self.exp_e = np.zeros(
            (int(
                self.Exp_Num_Nodes), int(
                self.Exp_Num_Nodes)))
        self.exp_M = np.zeros(
            (int(
                self.Exp_Num_Nodes), int(
                self.Exp_Num_Nodes)))
        self.exp_T = np.zeros(
            (int(
                self.Exp_Num_Nodes), int(
                self.Exp_Num_Nodes)))

        self.generate_neighborhood_matrix( PD_deck, PD_problem.exp_init_positions )
        print "Family:"
        print self.family
        self.compute_ext_state(PD_deck, PD_problem, t)
        print "Ext_state:"
        print self.exp_e
        self.compute_M(PD_problem, t)
        print "M:"
        print self.exp_M
        self.compute_T(PD_deck, PD_problem)
        self.compute_Ts(PD_deck)
        self.compute_exp_strain_energy(PD_problem, t)
        #self.compute_T(PD_deck, PD_problem, y)
        #self.compute_Ts(PD_deck, PD_problem)

    # Computes the deformations using the current positions and initial
    # positions of each node
    # Reminder: y[xi, t] = x[xi] + u[xi, t]
    def compute_ext_state(self, PD_deck, PD_problem, t):
        # Initialization for e
        for x_i in range(0, self.Exp_Num_Nodes):
            index_x_family = self.get_index_x_family(x_i)
            for x_p in index_x_family:
                self.exp_e[x_i, x_p] = np.absolute(
                    PD_problem.exp_displacement[t][x_p] - PD_problem.exp_displacement[t][x_i]) - np.absolute(PD_problem.exp_init_positions[
                        x_i] - PD_problem.exp_init_positions[x_p])

    def compute_M(self, PD_problem, t):
        M = np.zeros((int(self.Exp_Num_Nodes), int(self.Exp_Num_Nodes)))
        for x_i in range(0, self.Exp_Num_Nodes):
            y_i = (
                PD_problem.exp_displacement[t][x_i] +
                PD_problem.exp_init_positions[x_i])

            index_x_family = self.get_index_x_family(x_i)
            for x_p in index_x_family:
                y_p = (
                    PD_problem.exp_displacement[t][x_p] +
                    PD_problem.exp_init_positions[x_p])

                self.exp_M[x_i, x_p] = (y_p-y_i)/np.absolute(y_p-y_i)

    # Computes the weights for each PD node
    def weighted_function(self, PD_problem, PD_deck, x, x_i):
        Horizon = PD_problem.compute_horizon(PD_deck)
        index_x_family = self.get_index_x_family(x_i)
        Delta_V = PD_deck.Volume
        Influence_Function = PD_deck.Influence_Function
        result = 0
        for x_p in index_x_family:
            result = result + Influence_Function * \
                (PD_problem.exp_init_positions[x_p] - PD_problem.exp_init_positions[x_i]) ** 2 * Delta_V
        return result

    def compute_T(self, PD_deck, PD_problem):
        w = PD_deck.Influence_Function
        for x_i in range(0, self.Exp_Num_Nodes):
            for e_i in range(0, len(self.exp_e[x_i])):
                self.exp_T[x_i,
                           e_i] = self.exp_M[x_i,
                                             e_i] * (w / self.weighted_function(PD_problem, PD_deck,PD_problem.exp_init_positions[x_i], x_i)) * self.Modulus * self.exp_e[x_i, e_i]

    def compute_Ts(self, PD_deck):
        for x_i in range(0, self.Exp_Num_Nodes):
            if not(x_i - 1 < 0):
                self.exp_Ts[x_i] += self.exp_T[x_i, x_i - 1] - \
                    self.exp_T[x_i - 1, x_i]
            if not(x_i + 1 > self.Exp_Num_Nodes - 1):
                self.exp_Ts[x_i] += self.exp_T[x_i, x_i + 1] - \
                    self.exp_T[x_i + 1, x_i]

            #self.exp_Ts[x_i] = self.exp_Ts[x_i]* PD_deck.Volume
            self.exp_Ts[x_i] *= PD_deck.Volume

    def compute_exp_strain_energy(self, PD_problem, t):
        for x_i in range(0, self.Exp_Num_Nodes):
            self.exp_W[x_i] = self.exp_Ts[x_i] * \
                PD_problem.exp_displacement[t][x_i]

    def generate_neighborhood_matrix(self, PD_deck, x):
        self.family = np.zeros((int(
                self.Exp_Num_Nodes), int(
                self.Exp_Num_Nodes)))
        for x_i in range(0, len(x)):
            for x_p in range(0, len(x)):
                if x_p == x_i:
                    pass
                elif np.absolute(x_i - x_p) <= PD_deck.Horizon_Factor:
                    self.family[x_i][x_p] = 1
                else:
                    pass

    # Returns a list of addresses of the neighbors of a point x_i
    def get_index_x_family(self, x_i):
        #print (np.where( self.family[x_i] == 1 ))[0]
        return (np.where(self.family[x_i] == 1))[0]
