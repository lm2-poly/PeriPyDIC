# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca

from ..util import neighbor
from ..util import abstractions
import numpy as np

## Copmutes extension and force states out of displacement/nodes data obtained
# from the VIC3D CSV (DIC data)
class DIC_problem(abstractions.Problem):

    ## Constructor
    # Find neighbors for DIC data, computes the weighted function, then computes
    # actual positions, extension states and force states
    # @param deck Deck object containig data from the .yaml file
    def __init__(self, deck):
        ## NeighborSearch
        self.neighbors = neighbor.NeighborSearch(deck)

        # Compute the volume correction factor for each node
        self.compute_volume_correction(deck)

        # Compute the weighted volume for each node
        self.compute_weighted_volume(deck)

        ## Actual position from DIC result
        self.y = np.zeros((deck.num_nodes, deck.dim,2),dtype=np.float32)
        self.y[:,:,0] = deck.geometry.nodes[:,:]

        ## Internal forces
        self.force_int = np.zeros((deck.num_nodes, deck.dim,2),dtype=np.float32)

        ## Extension state
        self.ext = np.zeros( ( deck.num_nodes, self.neighbors.max_neighbors,2),dtype=np.float32 )
        
        if deck.critical_stretch > 0:
            ## Bond damage
            self.damage = np.zeros((deck.num_nodes,deck.num_nodes,deck.time_steps))


        if deck.material_type == "Elastic":
            from ..materials.elastic import Elastic_material
            mat_class = Elastic_material( deck, self, deck.geometry.act )
            self.update_force_data(mat_class)
            self.update_ext_state_data(mat_class)

        if not deck.damage_type == "None":    
            self.compute_damage(deck, self.y)
            self.update_pos(deck.geometry.act)
        
        self.update_damage(deck)    
        self.strain_energy = mat_class.strain_energy

    ## Records the force vector at each time step
    # @param mat_class Material class object for the elastic/viscoelastic material models
    def update_force_data(self, mat_class):
        ## Internal forces state
        self.force_int[:,:, 1] = mat_class.f_int


    ## Records the ext_state vector at each time step
    # @param mat_class Material class object for the elastic/viscoelastic material models
    def update_ext_state_data(self, mat_class):
        ## Extension state
        self.ext[:, :, 1] = mat_class.e

    ## Records the actual position vector at each time step
    # @param act Actual position obtained from DIC data
    def update_pos(self,act):
        ## Actual position state
        self.y[:,:, 1] = act
    
    ## Update the damage between bonds
    # @param     
    def update_damage(self):
        self.damage[:,:, 1] = self.neighbors.damage
