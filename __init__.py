# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 19:24:33 2015

@author: ilyass
"""

#from Class_1D-PD-elastic-problem import PD_1D_elastic_problem
from deck import PD_deck
from problem import PD_1D_problem
import numpy as np

logger = logging.getLogger(__name__)

data = PD_deck()

problem = PD_1D_problem( data )

#Build empty matrix u[row = node, column = time]
y = np.zeros( ( int(data.Num_Nodes), int(data.Num_TimeStep)) )
logger.info("Successfully created matrix u")
