# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 19:24:33 2015

@author: ilyass
"""

#from Class_1D-PD-elastic-problem import PD_1D_elastic_problem
from deck import PD_deck
from problem import PD_problem
from elastic import elastic_material
import numpy as np
import logging
import random


logger = logging.getLogger(__name__)

data = PD_deck()

problem = PD_problem( data )

y = np.zeros( ( int(data.Num_Nodes) ) )
for x_i in range(0, int(data.Num_Nodes)):
    y[x_i] = problem.x[x_i]
print 'Y init'
print y

print "X"
print problem.x

forces = elastic_material( data, problem, y )

print problem.quasi_static_solver( y, data, forces )

problem.write_data_to_csv(data, problem)

problem.plot_force(data)

force_plot = problem.plot_force(data)
force_plot.show()

position_plot = problem.plot_force(data)
position_plot.show()