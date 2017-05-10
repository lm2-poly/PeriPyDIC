# -*- coding: utf-8 -*-
#@author: ilyass.tabiai@polymtl.ca
#@author: rolland.delorme@polymtl.ca
#@author: patrick.diehl@polymtl.ca

from ..util import neighbor
from ..util import linalgebra
import numpy as np

## Copmutes extension and force states out of displacement/nodes data obtained
# from the VIC3D CSV (DIC data)
class DIC_problem():

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


        if deck.material_type == "Elastic":
            from ..materials.elastic import Elastic_material
            mat_class = Elastic_material( deck, self, deck.geometry.act )
            self.update_force_data(mat_class)
            self.update_ext_state_data(mat_class)
            

        self.update_pos(deck.geometry.act)
        self.strain_energy = mat_class.strain_energy

    ## Compute the volume correction factors
    # @param deck The input deck
    def compute_volume_correction(self,deck):
        ## Volume correction factor for each node
        self.volume_correction = np.ones( ( deck.num_nodes, self.neighbors.max_neighbors), dtype=np.float64 )
        for i in range(0, deck.num_nodes):
            index_x_family = self.neighbors.get_index_x_family(i)
            n = 0            
            for p in index_x_family:
                X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                r = deck.delta_X / 2.0
                if linalgebra.norm(X) > self.neighbors.horizon - r:
                    self.volume_correction[i,n] = (self.neighbors.horizon + r - linalgebra.norm(X)) / (deck.delta_X)
                else:
                    pass
                n += 1

    ## Compute the weighted volume for each node
    # @param deck The input deck
    def compute_weighted_volume(self, deck):
        ## Weighted volume for each node
        self.weighted_volume = np.zeros((deck.num_nodes),dtype=np.float64)
        for i in range(0, deck.num_nodes):
            index_x_family = self.neighbors.get_index_x_family(i)
            n = 0            
            for p in index_x_family:
                X = deck.geometry.nodes[p,:] - deck.geometry.nodes[i,:]
                self.weighted_volume[i] += deck.influence_function * (linalgebra.norm(X))**2 * self.volume_correction[i,n] * deck.geometry.volumes[p]
                n += 1

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
