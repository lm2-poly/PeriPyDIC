# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 16:16:07 2015

@author: ilyass.tabiai@polymtl.ca
@author: rolland.delorme@polymtl.ca
@author: patrick.diehl@polymtl.ca
"""

import numpy as np

class viscoelastic_material():
    
    def __init__(self, PD_deck, PD_problem, y, t_n):
        self.len_x = PD_deck.num_nodes_x
        self.Relax_Modulus = PD_deck.relax_modulus
        self.Relax_Time = PD_deck.relax_time
        self.w = PD_deck.influence_function
        self.M = PD_problem.compute_m(y)
        self.compute_ext_state(PD_deck, PD_problem, y)
        self.compute_ext_state_visco(PD_deck, PD_problem, y, t_n )   
        self.compute_t_visco(PD_deck, PD_problem, y)
        self.compute_T(PD_deck, PD_problem, y)
        self.compute_Ts(PD_deck, PD_problem)
               
    #Computes the deformations using the current positions and initial 
    #positions of each node
    #Reminder: y[xi, t] = x[xi] + u[xi, t]
    def compute_ext_state(self, PD_deck, PD_problem, y):        
        e = np.zeros( (self.len_x, self.len_x) )
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family(x_i)
            for x_p in index_x_family: 
                e[x_i, x_p] = np.absolute(y[x_p] - y[x_i]) - np.absolute(PD_deck.geometry.pos_x[x_p] - PD_deck.geometry.pos_x[x_i])       
        self.e = e

    def compute_ext_state_visco(self, PD_deck, PD_problem, y, t_n):        

        e_visco = np.zeros( (self.len_x, self.len_x, len(self.Relax_Time)) )
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family(x_i)
            for x_p in index_x_family: 
                for k in range(1, len(self.Relax_Time)):
                    temp_exp = np.exp((- PD_deck.delta_t)/(self.Relax_Time[k]))
                    delta_e = self.e[x_i, x_p] - PD_problem.ext[x_i, x_p, t_n-1]
                    beta = 1-(self.Relax_Time[k]*(1-temp_exp)) / PD_deck.delta_t
                    e_visco[x_i, x_p, k] = PD_problem.ext[x_i, x_p, t_n-1] *(1-temp_exp) + PD_problem.ext_visco[x_i, x_p, k, t_n-1]*temp_exp + beta*delta_e         
        self.e_visco = e_visco

    def compute_t_visco(self, PD_deck, PD_problem, y):
        t_visco = np.zeros( (self.len_x, self.len_x ) )
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family(x_i)
            for x_p in index_x_family: 
                for k in range(1, len(self.Relax_Time)):
                    t_visco[x_i, x_p] = t_visco[x_i, x_p] + (self.w / PD_problem.weighted_function(PD_deck, PD_deck.geometry.pos_x, x_i))*self.Relax_Modulus[k]*(self.e[x_i, x_p] - self.e_visco[x_i, x_p, k])
        self.t_visco = t_visco

    def compute_T(self, PD_deck, PD_problem, y):      
        tscal = np.zeros( (self.len_x, self.len_x ) )
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family(x_i)
            for x_p in index_x_family:
                tscal[x_i, x_p] = (self.w / PD_problem.weighted_function(PD_deck, PD_deck.geometry.pos_x, x_i))*self.Relax_Modulus[0]*self.e[x_i, x_p] + self.t_visco[x_i, x_p]
        self.tscal = tscal
        
        T = np.zeros( (self.len_x, self.len_x ) )
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family(x_i)
            for x_p in index_x_family:
                T[x_i, x_p] = tscal[x_i, x_p] * self.M[x_i, x_p]
        self.T = T
        
    def compute_Ts(self, PD_deck, PD_problem):
        Ts = np.zeros( (self.len_x ) )
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family(x_i)
            for x_p in index_x_family:
                Ts[x_i] = Ts[x_i] + self.T[x_i, x_p] - self.T[x_p, x_i]
            Ts[x_i] = Ts[x_i] * PD_deck.geometry.volumes[x_i]
        self.Ts = Ts