# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 16:16:07 2015

@author: ilyass.tabiai@gmail.com
@author: rolland.delorme@gmail.com
"""

from deck_visco import PD_deck
#from problem import PD_problem
from problem_sym import PD_problem
import numpy as np

class viscoelastic_material():
    
    def __init__(self, PD_deck, PD_problem, y, t_n):
        self.len_x = len(PD_problem.x)
        self.Modulus, self.Relaxation_Time = PD_deck.get_viscoelastic_material_properties()
        self.compute_ext_state(PD_deck, PD_problem, y)
        self.compute_ext_state_visco(PD_deck, PD_problem, y, t_n )   
        self.compute_t_visco(PD_deck, PD_problem, y)
        self.compute_T(PD_deck, PD_problem, y)
        self.compute_Ts(PD_deck, PD_problem)
        
        
    #Computes the deformations using the current positions and initial 
    #positions of each node
    #Reminder: y[xi, t] = x[xi] + u[xi, t]
    def compute_ext_state(self, PD_deck, PD_problem, y):        
        e = np.zeros( (int(PD_deck.Num_Nodes), int(PD_deck.Num_Nodes)) )

        for x_i in range(0, len(PD_problem.x)):
            index_x_family = PD_problem.get_index_x_family(PD_problem.x, x_i)
            #print index_x_family
            for x_p in index_x_family: 
                e[x_i, x_p] = np.absolute(y[x_p] - y[x_i]) - np.absolute(PD_problem.x[x_p] - PD_problem.x[x_i])       
        self.e = e

    def compute_ext_state_visco(self, PD_deck, PD_problem, y, t_n):        

        e_visco = np.zeros( (int(PD_deck.Num_Nodes), int(PD_deck.Num_Nodes), len(self.Relaxation_Time)) )
        for x_i in range(0, len(PD_problem.x)):
            index_x_family = PD_problem.get_index_x_family(PD_problem.x, x_i)
            for x_p in index_x_family: 
                for k in range(1, len(self.Relaxation_Time)):
                    temp_exp = np.exp((- PD_deck.Delta_t)/(self.Relaxation_Time[k]))
                    delta_e = self.e[x_i, x_p] - PD_problem.ext[x_i, x_p, t_n-1]
                    beta = 1-(self.Relaxation_Time[k]*(1-temp_exp)) / PD_deck.Delta_t
                    e_visco[x_i, x_p, k] = PD_problem.ext[x_i, x_p, t_n-1] *(1-temp_exp) + PD_problem.ext_visco[x_i, x_p, k, t_n-1]*temp_exp + beta*delta_e         
        self.e_visco = e_visco

    def compute_t_visco(self, PD_deck, PD_problem, y):
        w = PD_deck.Influence_Function
        M = PD_problem.compute_m(PD_deck.Num_Nodes, y)
        
        t_visco = np.zeros( (int(PD_deck.Num_Nodes), int(PD_deck.Num_Nodes) ) )
        for x_i in range(0, len(PD_problem.x)):
            index_x_family = PD_problem.get_index_x_family(PD_problem.x, x_i)
            for x_p in index_x_family: 
                for k in range(1, len(self.Relaxation_Time)):
                    t_visco[x_i, x_p] = t_visco[x_i, x_p] + (w/PD_problem.weighted_function(PD_deck, PD_problem.x, x_i))*self.Modulus[k]*(self.e[x_i, x_p] - self.e_visco[x_i, x_p, k])
        self.t_visco = t_visco

    def compute_T(self, PD_deck, PD_problem, y):
        w = PD_deck.Influence_Function
        M = PD_problem.compute_m(PD_deck.Num_Nodes, y)
        
        tscal = np.zeros( (int(PD_deck.Num_Nodes), int(PD_deck.Num_Nodes) ) )
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family( PD_problem.x, x_i)
            for x_p in index_x_family:
                tscal[x_i, x_p] = (w/PD_problem.weighted_function(PD_deck, PD_problem.x, x_i))*self.Modulus[0]*self.e[x_i, x_p] + self.t_visco[x_i, x_p]
        self.tscal = tscal
        
        T = np.zeros( (int(PD_deck.Num_Nodes), int(PD_deck.Num_Nodes) ) )
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family( PD_problem.x, x_i)
            for x_p in index_x_family:
                T[x_i, x_p] = tscal[x_i, x_p] * M[x_i, x_p]
        self.T = T
        
    def compute_Ts(self, PD_deck, PD_problem):
        Ts = np.zeros( (int(PD_deck.Num_Nodes) ) )
        for x_i in range(0, self.len_x):
            index_x_family = PD_problem.get_index_x_family( PD_problem.x, x_i)
            for x_p in index_x_family:
                Ts[x_i] = Ts[x_i] + self.T[x_i, x_p] - self.T[x_p, x_i]
            Ts[x_i] = Ts[x_i] * PD_deck.Volume
        self.Ts = Ts

