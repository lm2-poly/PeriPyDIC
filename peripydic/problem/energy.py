#-*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca

import numpy as np
from scipy import linalg

from ..util import neighbor
from ..util import linalgebra
from ..util import abstractions

import sys
from peripydic.IO import deck

class Energy_problem(abstractions.Problem):
    
    ## Constructor
    # @param deck The input deck
    def __init__(self, deck):
        
        ## NeighborSearch
        self.neighbors = neighbor.NeighborSearch(deck)

        # Compute the volume correction factor for each node
        self.compute_volume_correction(deck)

        # Compute the weighted volume for each node
        self.compute_weighted_volume(deck)
        
        
    def jacobian_matrix(self, deck, y, p):
        eps = deck.solver_perturbation
        jacobian = np.zeros((deck.compare_length , deck.dim),dtype=np.float64)
        #print p
        index = 0
        
        for i in deck.nodes_compare:
            for j in range(0,len(p.shape)):
                
                if deck.material_type == "Elastic":
                    if deck.dim == 1:
                        from ..materials.elastic import Elastic_material
                        deck.young_modulus = p[j] + eps
                        mat_class_p = Elastic_material( deck, self, y )

                        deck.young_modulus = p[j] - eps
                        mat_class_m = Elastic_material( deck, self, y )
                        jacobian[index][j] = (mat_class_p.strain_energy[i] - mat_class_m.strain_energy[i]) / (2. * eps)
                        
                        index +=1
        return jacobian
    
    def newton_step(self, deck, y,p):
        
        jacobian = self.jacobian_matrix(deck, y,p)
      
        S = linalg.pinv(jacobian)
        
        energy = np.zeros((deck.compare_length,deck.dim))
        
        for i in range(0,len(energy)):
            for j in range(0,deck.dim):
                energy[i][j] = deck.measured_energy - jacobian[i][j]
        
     
        if deck.dim == 1:
            p[0] =  np.dot(S,energy)[0]
        if deck.dim == 2:
            res = np.dot(S,energy)
            p[0] = res[0]
            p[1] = res[1]
       
        for i in range(0,len(energy)):
            from ..materials.elastic import Elastic_material
            if deck.dim == 1:
                deck.young_modulus = p[0]
            if deck.dim == 2:
                deck.bulk_modulus = p[0]
                deck.shear_modulus = p[1]
                
            mat_class = Elastic_material( deck, self, y )
            energy[i] = deck.measured_energy - mat_class.strain_energy[deck.nodes_compare[i]]
           
        return linalgebra.norm(energy) 
    
    def solver(self,deck):
        
        
        if deck.dim == 1:
            p = np.zeros((1))
            p.fill(deck.young_modulus)
            
        if deck.dim == 2:
            p = np.zeros((2))
            p[0] = deck.bulk_modulus
            p[1] = deck.shear_modulus
            
        
        res = float('inf')
        iteration = 1
      
        while res >= deck.solver_tolerance and iteration <= deck.solver_max_it :
            
            res = self.newton_step(deck, deck.geometry.act,p)
            print iteration , res 
            iteration += 1
            
        print p
        
        