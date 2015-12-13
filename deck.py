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
                self.initial_data = initial_data['data']
                self.read_data(initial_data['data'])
                self.compute_volumes()
            except:
                logger.error("The XML file is broken")
                
    def read_data(self, initial_data):
        self.Horizon_Factor = int(initial_data['Discretization']['Horizon_Factor'])
        self.N_Delta_t = int(initial_data['Discretization']['N_Delta_t'])
        self.Final_Time = float(initial_data['Discretization']['Final_Time'])
        self.Delta_t = float(self.Final_Time/self.N_Delta_t)
        self.N_Delta_x = int(initial_data['Discretization']['N_Delta_x'])
        #Total length of the 1D bar
        self.Length = float(initial_data['Geometry']['Length'])
        #Number of PD nodes "meshing" the bar
        self.Num_Nodes = int(self.N_Delta_x + 1 + 2 * self.Horizon_Factor)
        #Distance between each couple of PD nodes
        self.Delta_x = float(self.Length / self.N_Delta_x)
        self.Length_Tot = self.Length + 2 * self.Horizon_Factor * self.Delta_x
        #Compute the total number of timesteps
        self.Num_TimeStep = int(self.N_Delta_t + 1)
        self.Loading_Flag = str(initial_data['Boundary_Conditions']['Type'])
        
    
    def compute_volumes(self):
        self.Surface = float(self.initial_data['Geometry']['Surface'])
        self.Volume = self.Surface * self.Delta_x
        #Compute the force density on an elementary volume
        self.Volume_Boundary = self.Volume * self.Horizon_Factor
          
    def compute_force_density(self):
        self.Force_Density = self.Force/self.Volume_Boundary

    def get_parameters_loading_ramp(self):
        self.Ramp_Time = float(self.initial_data['Boundary_Conditions']['Ramp_Time'])
        self.Force = float(self.initial_data['Boundary_Conditions']['Force'])
        self.compute_force_density()
        