# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 18:11:02 2015

@author: ilyass.tabiai@gmail.com
"""
import soft_utilities
import funct_utilities
import logging
from scipy.optimize import fsolve
import timeit
from Class_PD_deck import PD_deck

logger = logging.getLogger(__name__)

class PD_1D_elastic_problem():
    
    def __init__(self):
        #Import initial data
        initial_data = PD_deck()
