# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 18:49:10 2015

@author: ilyass
"""

import xmltodict
import logging
from xml.parsers.expat import ExpatError
import numpy as np

logger = logging.getLogger(__name__)

class PD_deck():
    
    def __init__(self):
        with open("deck.xml") as deck:
            try:
                initial_data = xmltodict.parse( deck.read() )
                self.read_data(initial_data['data'])
                self.get_pd_nodes(initial_data['data'])
            except:
                logger.error("The XML file is broken")
                
    def read_data(self, initial_data):
        self.Horizon_Factor = int(initial_data['Discretization']['Horizon_Factor'])
        self.N_Delta_t = int(initial_data['Discretization']['N_Delta_t'])
        self.Final_Time = float(initial_data['Discretization']['Final_Time'])
        self.N_Delta_x = int(initial_data['Discretization']['N_Delta_x'])
        self.Length = float(initial_data['Geometry']['Length'])
        self.Num_Nodes = int(self.N_Delta_x + 1 + 2 * self.Horizon_Factor)
        self.Delta_x = float(self.Length / self.N_Delta_x)
        self.Length_Tot = self.Length + 2 * self.Horizon_Factor * self.Delta_x
        self.num_TimeStep = int(self.N_Delta_t + 1)
    def get_pd_nodes(self, initial_data):
        # Define x
        x = np.zeros( self.Num_Nodes )    
        for i in range(0, self.Num_Nodes):
            x[i] = -self.Length_Tot/2 + i * self.Delta_x
        self.x = x