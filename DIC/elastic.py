#-*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca

#from problem.pd import PD_problem
import numpy as np

class elastic_material_dic():

    def __init__(self, deck, problem):
        ## Scalar influence function
        self.w = deck.influence_function
        
        ## Weighted volume
        self.Weighted_Volume = problem.weighted_volume

        if deck.dim == 1:
            ## Young modulus of the material
            self.Young_Modulus = deck.young_modulus
        if deck.dim == 2:
            ## Bulk modulus of the material
            self.K = deck.bulk_modulus
            ## Shear modulus of the material
            self.Mu = deck.shear_modulus
            ## Poisson ratio of the material
            self.Nu = (3. * self.K - 2. * self.Mu) / (2. * (3. * self.K + self.Mu))

        if deck.dim == 3:
            ## Bulk modulus of the material
            self.K = deck.bulk_modulus
            ## Shear modulus of the material
            self.Mu = deck.shear_modulus
