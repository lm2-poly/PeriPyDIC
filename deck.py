# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 18:49:10 2015

@author: ilyass.tabiai@gmail.com
@author: rolland.delorme@gmail.com
"""

import xmltodict
import logging
from xml.parsers.expat import ExpatError
import numpy as np

logger = logging.getLogger(__name__)

#This class is pretty self-explanatory, it retrieves data from the XML deck and
#records it as variables of the PD_deck class object
class PD_deck():
    
    def __init__(self, path):
        with open(path) as deck:
            try:
                #Here we use the method parse from xmltodict package
                initial_data = xmltodict.parse( deck.read() )
                self.initial_data = initial_data['data']
                self.read_data(initial_data['data'])
                self.compute_volumes()
            except:
                logger.error("The XML file is broken")
                
    def read_data(self, initial_data):
        self.Horizon_Factor = int(initial_data['Discretization']['Horizon_Factor'])
        self.N_Steps_t = int(initial_data['Discretization']['N_Steps_t'])
        self.Final_Time = float(initial_data['Discretization']['Final_Time'])
        self.Delta_t = float(self.Final_Time/self.N_Steps_t)
        self.N_Nodes_Bar = int(initial_data['Discretization']['N_Nodes_Bar'])
        #Length of the 1D bar
        self.Length_Bar = float(initial_data['Geometry']['Length_Bar'])
        #Number of PD nodes "meshing" the bar including the boundary nodes
        self.Num_Nodes = int(self.N_Nodes_Bar + 2 * self.Horizon_Factor)
        #Distance between each couple of PD nodes
        self.Delta_x = float(self.Length_Bar / self.N_Nodes_Bar)
        #Lenght of the model including the boundary nodes        
        self.Length_Tot = self.Num_Nodes * self.Delta_x
        #Compute the total number of timesteps
        self.Num_TimeStep = int(self.N_Steps_t + 1)
        self.Loading_Flag = str(initial_data['Boundary_Conditions']['Type'])
        self.Material_Flag = str(initial_data['Material']['Type'])
        self.Influence_Function = float(initial_data['Discretization']['Influence_Function'])
    
    #We compute the volumes in order to be able to compute stresses later on
    #To make comparisons with continuum mechanics easier
    def compute_volumes(self):
        self.Surface = float(self.initial_data['Geometry']['Surface_Bar'])
        self.Volume = self.Surface * self.Delta_x
        self.Volume_Boundary = self.Volume * self.Horizon_Factor
        
    #Only called if the Type of Boundary Conditions in the XML deck is RAMP
    def get_parameters_loading_ramp(self):
        self.Ramp_Time = float(self.initial_data['Boundary_Conditions']['Ramp_Time'])
        self.Force = float(self.initial_data['Boundary_Conditions']['Force'])
        self.compute_force_density()

    #Compute the force density on an elementary volume
    def compute_force_density(self):
        self.Force_Density = self.Force/self.Volume_Boundary
        
    #Only called it the Type of material in the XML deck is ELASTIC
    def get_elastic_material_properties(self):
        Modulus = float(self.initial_data['Material']['E_Modulus'])  
        return Modulus 