# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 19:24:33 2015

@author: ilyass.tabiai@gmail.com
@author: rolland.delorme@gmail.com
"""
#from Class_1D-PD-elastic-problem import PD_1D_elastic_problem
import logging
logger = logging.getLogger(__name__)

import numpy as np
import random

#Load the PD_deck class and create a PD_deck object
#from deck_elas import PD_deck
from deck_visco import PD_deck
data = PD_deck()

#Load the PD_problem class and create a PD_problem object
from problem import PD_problem
problem = PD_problem( data )

#Create an initial guess vector here based on a linear distribtuion of the 
#nodes, disturbed by a small random coefficient using a method provided in the
#PD_problem class
x_0 = problem.random_initial_guess( problem.x, data ) 
#x_0 is our initial guess
#print x_0

#Load the elastic_material class and compute first step PD forces
#from elastic import elastic_material
#forces = elastic_material( data, problem, x_0 )

#Solve the problem
problem.quasi_static_solver( x_0, data )
problem.strain_center_bar( data )

#Check the position of PD nodes at each time step
print np.around(problem.y,decimals=5)

#Check the PD force value at each node at the 5th time step
#print problem.forces[:, 5]

#Check the strain in the middle of the bar at each time step
print np.around(problem.strain,decimals=10)
#Write the results to a CSV file
problem.write_data_to_csv(data)
#The problem resolution (time step by time step) is now written in a 
#csv file called data_csv in the current folder

