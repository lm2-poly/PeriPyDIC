# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca

import util.neighbor
from DIC import elastic
import numpy as np

class DIC_problem():
    
    def __init__(self, deck):
        ## Amount of nodes
        #self.len_x = deck.dim *  deck.num_nodes
        ## Neighbors for each DIC pixel are stored here
        #self.neighbors = util.neighbor.NeighborSearch(deck)
        
        ## Actual positions
        #self.y = np.zeros( ( self.len_x, 2) )
        
        ## Force 
        #self.force_int = np.zeros((deck.num_nodes, deck.dim),dtype=np.float64)
        
        #self.material = elastic.elastic_material_dic(deck,self)