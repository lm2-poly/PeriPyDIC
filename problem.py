# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 12:50:28 2015

@author: ilyass
"""

import logging
from scipy.optimize import fsolve
import timeit
from Class_PD_deck import PD_deck
import numpy as np

logger = logging.getLogger(__name__)

class PD_1D_problem():
    
    def __init__(self, PD_deck):
        #Import initial data
        self.get_pd_nodes(PD_deck)
        self.compute_b(PD_deck)
        
        
    def compute_b(self, PD_deck):       
        #Build  matrix b[row = node, column = time]
        b = np.zeros( ( int(PD_deck.Num_Nodes), int(PD_deck.Num_TimeStep)) )   
        PD_deck.get_parameters_loading_ramp()
        if PD_deck.Loading_Flag == "RAMP":
            PD_deck.get_parameters_loading_ramp()
            for x_i in range(0, PD_deck.Horizon_Factor):
                for t_n in range(1, int(PD_deck.Num_TimeStep)):
                    b[x_i, t_n] = - self.ramp_loading( PD_deck, t_n )
            for x_i in range(len(self.x) - PD_deck.Horizon_Factor, len(self.x) ):
                for t_n in range(1, int(PD_deck.Num_TimeStep)):
                    b[x_i, t_n] = self.ramp_loading( PD_deck, t_n )
        else:
            logger.error("There is a problem with the Boundary Conditions in your XML deck.")
        self.b = b
        
    def get_pd_nodes(self, PD_deck):
        # Define x
        x = np.zeros( PD_deck.Num_Nodes )    
        for i in range(0, PD_deck.Num_Nodes):
            x[i] = -PD_deck.Length_Tot/2 + i * PD_deck.Delta_x
        self.x = x
        
    def ramp_loading(self, PD_deck, t_n):     
        Time_t = PD_deck.Delta_t*t_n
        for x_i in range(0, int(PD_deck.Num_Nodes)):
            if Time_t <= PD_deck.Ramp_Time:
                result = (PD_deck.Force_Density*Time_t)/PD_deck.Ramp_Time   
                return result
            else:
                result = PD_deck.Force_Density
                return result